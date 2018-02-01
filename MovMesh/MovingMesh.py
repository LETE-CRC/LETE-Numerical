#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 2015 by helio
Modified on Tue Nov 01 2017 by filipi
"""

import numpy as np
import re
from scipy.sparse.linalg import bicg # bicgstab
from timeit import default_timer as timer
import os
import sys

np.set_printoptions(threshold=sys.maxint)

#==============================================================================
# --------------------- INPUTS E SELEÇÕES NECESSÁRIAS ------------------------#
#==============================================================================
LaplacianScheme = 'linear'
gamma = 1000000 # RIGIDEZ DA MALHA
print(gamma)
CI='deformation125' #condições iniciais
scalar = 'quality'
tinitial = 0
dt = 1
tfinal = 2
###############################################################################
# --------------- ---------------- FUNÇÕES ---------------------------------- #
###############################################################################
    
def read_file(arquivo): # funcao que le arquivos de malha do OpenFOAM
    quantidade = 0
    with open(arquivo, 'r') as infile: # abre o arquivo a ser lido
    
        for line in infile:
            if line.strip() == "(": # quando encontra (,
                break  # remove tudo anterior ao primeiro (
            quantidade = line.strip() # guarda o valor (str) da qtd de pontos, faces, etc

        quantidade = int(quantidade)# qtd em forma de numero
        arquivo = [] # declara lista de pontos, faces,etc
        
        for line in infile:
            if line.strip() == ")": # quando encontra ) sozinho na linha, para(break)
                break
            line = line.rsplit('(', 1)[-1];line = line.strip()
            line = line.translate(None, ")");line = line.split(' ') # formatacoes nas linhas
            arquivo.append(line) # adiciona a linha na lista
            
        return (arquivo, quantidade)

#=============================================================================#

def read_scalarlist(path): # func que le arquivos de escalares (owner/neighbour)
    with open(path, 'r') as infile:
        arquivo = infile.read()
        arquivo = arquivo[arquivo.find("(")+1:arquivo.find(")")]
        lista = re.findall(r'\s*(\d*\.\d+|\d+)', arquivo)
    return lista

#=============================================================================#

def read_bcs_mesh(): # funcao que le condicoes de contorno da malha do OpenFOAM
    with open(dirMalhaBoundary, 'r') as infile: # abre o arquivo a ser lido
        arquivo = infile.read() # le o arquivo como unica string
        arquivo = ''.join(arquivo.split()) # retira todos espacos em branco
        namebcs = re.findall(r'(\w+)\{type.\}*', arquivo) # encontra nomes das bcs
        tipobc = re.findall(r'type(\w+)', arquivo) # encontra tipo
        nfaces = re.findall(r'nFaces(\w+)', arquivo) # encontra id primeira face na lista de faces
        startface = re.findall(r'startFace(\w+)', arquivo) # encontra id primeira face na lista de faces
        nfaces = np.asarray(nfaces);nfaces = nfaces.astype(np.int)
    return (namebcs, tipobc, nfaces, startface)

#=============================================================================#

def write_OF(dirFinal,varName,solution): # func que escreve resultados

    pathFinal = dirFinal + varName
    if not os.path.exists(dirFinal):
        os.makedirs(dirFinal)
    
    res = str(solution).replace(']', ')')
    #correção de espaços, parentesis e pontos    
    res = res.replace('[','(')
    res = res.replace('((','(')
    res = res.replace('))',')')
    res = res.replace(' (','(')
    res = res.replace('( ','(')
    res = res.replace('.)',')')
    res = res.replace('. ',' ')
    res = res.replace('  ',' ')
    res = res.replace('  ',' ')
    res = res.replace('  ',' ')
    res = ''.join(res)
    
    header = '/*--------------------------------*- C++ -*----------------------------------*\\'
    header += '\n| =========                 |                                                 |'
    header += '\n| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'
    header += '\n|  \\\    /   O peration     | Version:  5.0                                   |'
    header += '\n|   \\\  /    A nd           | Web:      www.OpenFOAM.org                      |'
    header += '\n|    \\\/     M anipulation  |                                                 |'
    header += '\n\*---------------------------------------------------------------------------*/'
    header += '\nFoamFile'
    header += '\n{'
    header += '\n    version     2.0;'
    header += '\n    format      ascii;'
    header += '\n    class       vectorField;'
    header += '\n    location    "%s";' % dirFinalr
    header += '\n    object      %s;' % varName
    header += '\n}'
    header += '\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
    header += '\n'
    header += '\n'
    header += '\n%d' % Npoints
    header += '\n(\n'
    
    end = '\n)'
    end += '\n'
    end += '\n'
    end += '\n// ************************************************************************* //'

    with open(pathFinal, 'w') as outfile:
            outfile.write(header)
            for item in res:            
                outfile.write(item)
            outfile.write(end)

    #arrumando outputs parte 1
    with open(pathFinal) as f:
    	newText=f.read().replace('  ', ' ')
	
    with open(pathFinal, "w") as f:
  	  f.write(newText)

    #arrumando outputs parte 2
    with open(pathFinal) as f:
    	newText=f.read().replace(' )', ')')
	
    with open(pathFinal, "w") as f:
  	  f.write(newText)
    
    return ()
#=============================================================================#

def calc_mesh_faces():
    areaFace = np.zeros(len(faces))
    areaFaceV = np.zeros((len(faces),3))
    cFace = np.zeros((len(faces),3)) # center of faces
    
    for i in range(len(faces)): # loop em todas as faces
        
        for j in range(2, len(faces[i])): #varre os indices de cada face (i) comecando pelo terceiro

            vetorA = points[faces[i][j-1]] - points[faces[i][j-2]] #os dois vetores para o prodVetorial
            vetorB = points[faces[i][j]] - points[faces[i][j-2]] #...
            # vetor area atraves da soma de cada vet area triangulo
            areaFaceV[i] = areaFaceV[i] + np.cross(vetorA, vetorB)/2 
            
        areaFace[i] = np.linalg.norm(areaFaceV[i]) # area de cada face em modulo
        
        # soma e divisao pela qtd de pontos para centro de cada face (i)
        cFace[i] = (np.sum(points[faces[i]],axis=0))/len(faces[i])
        
    return (areaFace, areaFaceV, cFace)

#=============================================================================#

def calc_mesh_vol():
    Nvolumes = np.max(owner)+1 #+1 contagem inicia em 0...
    cVol = np.zeros((Nvolumes,3)); vol = np.zeros(Nvolumes)
    qtdFacesVol = np.zeros(Nvolumes)
    
    for i in range(Nvolumes): # loop em todos os volumes
        # identifica as faces que pertencem ao volume (i)
        facesofvolumeO = np.where(owner==i)[0]
        # idem so que na lista neighbour pois face nomeada uma unica vez nas listas
        facesofvolumeN = np.where(neighbour==i)[0]
        # quantidade de faces em cada volume
        qtdFacesVol[i] = qtdFacesVol[i] + len(facesofvolumeO) + len(facesofvolumeN)
        
        # centro volume soma nos das faces owner e neighbour / total nos
        cVol[i] = (np.sum(cFace[facesofvolumeO],0) + np.sum(cFace[facesofvolumeN],0))/qtdFacesVol[i]

        for j in facesofvolumeO:
            vol[i] += np.abs(np.dot(cVol[i]-cFace[j], areaFaceV[j]))
    
        for j in facesofvolumeN:
            vol[i] += np.abs(np.dot(cVol[i]-cFace[j], areaFaceV[j]))
    
        vol[i] = vol[i]/3
    return (Nvolumes, cVol, vol)

#=============================================================================#

def solucao():
    A = np.zeros([Nvolumes,Nvolumes]);Su = np.zeros(Nvolumes) ## termo FONTE
    Sp = np.zeros(Nvolumes); SuBC = np.zeros(Nfaces); SpBC = np.zeros(Nfaces)
    d = np.zeros(Nneighbour);dVet = np.zeros((Nneighbour,3))
    dPfV = np.zeros((Nfaces,3)); dPf = np.zeros((Nfaces))
    dDeltaV = np.zeros((Nneighbour,3));dDelta = np.zeros(Nneighbour)
    aNf = np.zeros(Nneighbour)

    for i in range(Nneighbour): # loop faces internas
        o = owner[i] ; n = neighbour[i]
        dVet[i] = cVol[n] - cVol[o] # distancia entre centros vol P e N
        d[i] = np.linalg.norm(dVet[i]) # modulo distancia entre centros
        dDeltaV[i] = np.cross(dVet[i], areaFaceV[i]) # corecao nao-ortogonalidade

        if np.linalg.norm(dDeltaV[i])<=1e-5: # caso malha ortogonal dDeltaV muito pequeno
            dDeltaV[i] = areaFaceV[i] # e portanto dDeltaV e igual ao vetor area
        else:
            dDeltaV[i] = dVet[i]*np.square(areaFace[i])/(np.dot(dVet[i],areaFaceV[i])) # se nao, aplicar correcao
        
        dDelta[i] = np.linalg.norm(dDeltaV[i]) # modulo dist corrigida ortogonalidade
        
        ## Definicao dos esquemas de discretizacao do laplaciano
        if LaplacianScheme == 'linear':
            # coef face interna ******** discretizacao linear ************
            Laplacian = dDelta[i] * gamma / d[i]
            
        elif LaplacianScheme =='none':
            Laplacian = 0

        aNf[i] = Laplacian

        A[(o,n)] += - aNf[i] # posicionamento aNf superior
        A[(n,o)] += - aNf[i] # idem inferior

    for i in range(Nfaces):
        o = owner[i]
        dPfV[i] = cFace[i] - cVol[o] # distancia do centro volume ate o centro da face
        dPf[i] = np.linalg.norm(dPfV[i])

    with open(dirInicial, 'r') as infile:
        arquivo = infile.read()
        
    for i in NameBCs:
        my_regex1 = re.escape(i) + r"\s*\{\s*type\s*(\w+)\s*\}*"
        my_regex2 = re.escape(i) + r"\s*\{\s*type\s*\w+\;\s*value\s*uniform\s*\((.+)\)"
        TipobcDic = re.findall(my_regex1, arquivo)
        TipobcDic = str(TipobcDic).strip('[\']')
        valor = re.findall(my_regex2, arquivo, re.MULTILINE)#le o arquivo 
        valor = str(valor).strip('[\']')                    #converte para string  
        valor=np.fromstring(valor, dtype=float, sep=" ")    #converte para float
        try:
            valor = float(valor)
        except:
            pass
            
        #print re.escape(i),'--->', TipobcDic, '--->',valor # Mostra no terminal as condicoes de contorno
        
        if TipobcDic == 'empty' or TipobcDic == 'symmetryPlane':
            j = NameBCs.index(i)
            for k in range(int(startFace[j]), int(startFace[j])+int(nFaces[j])):
                SuBC[k] = 0
                SpBC[k] = 0
            
        if TipobcDic == 'fixedValue':
            j = NameBCs.index(i)
            for k in range(int(startFace[j]), int(startFace[j])+int(nFaces[j])):
                SuBC[k] = (areaFace[k] * gamma / dPf[k]) * valor[dim]
                SpBC[k] = -areaFace[k] * gamma / dPf[k]

        
    for i in range(Nfaces):
        o = owner[i] # para cada face, o é o valor da lista owner (n do vol de controle)
        # Su (lado direito) e o inicial (com as geracoes volumetricas) mais contribuicoes cada face
        Su[o] = Su[o] + SuBC[i]
        Sp[o] = Sp[o] + SpBC[i]

    for i in range(len(A)):
        A[i,i] = -np.sum(A[i])-Sp[i] # aP da soma dos vizinhos
        
    return (A, Su)

#=============================================================================#

def interpolator_corrector():
    dCent = np.zeros(Nneighbour)
    dVolN = np.zeros(Nneighbour)
    dVolO = np.zeros(Nneighbour)
    OutputFaces = np.zeros([Nowner,3])
    DeslocamentoA = np.zeros([Npoints,3])
    DeslocamentoV = np.zeros([Npoints,3])
    DeslocamentoBC = np.empty([Nfaces,3])
    DeslocamentoBC[:]=np.nan
    
    for i in range(Nneighbour): # loop faces internas
            o = owner[i] ; n = neighbour[i]
            dCent[(i)] = np.linalg.norm(cVol[n] - cVol[o]) # módulo da distância entre centros vol N e O
            dVolN[(i)] = np.linalg.norm(cVol[n] - cFace[i] ) # módulo da distância entre centros vol N e face
            dVolO[(i)] = np.linalg.norm(cFace[i] - cVol[o]) # módulo da distância entre centros vol face e O
            OutputFaces[i,dim] = (Output[o]*dVolO[i]+Output[n]*dVolN[i])/dCent[i] #interpolando o resultado para as faces
    
    with open(dirInicial, 'r') as infile:
        arquivo = infile.read()
        
    for i in NameBCs:
            my_regex1 = re.escape(i) + r"\s*\{\s*type\s*(\w+)\s*\}*"
            my_regex2 = re.escape(i) + r"\s*\{\s*type\s*\w+\;\s*value\s*uniform\s*\((.+)\)"
            TipobcDic = re.findall(my_regex1, arquivo)
            TipobcDic = str(TipobcDic).strip('[\']')
            valor = re.findall(my_regex2, arquivo, re.MULTILINE)#le o arquivo 
            valor = str(valor).strip('[\']')                    #converte para string  
            valor=np.fromstring(valor, dtype=float, sep=" ")    #converte para float
            try:
                valor = float(valor)
            except:
                pass
                
            if TipobcDic == 'empty' or TipobcDic == 'symmetryPlane':
                j = NameBCs.index(i)
                for k in range(int(startFace[j]), int(startFace[j])+int(nFaces[j])):
                    o = owner[k]
                    OutputFaces[k,dim] = Output[o]             
                    
            if TipobcDic == 'fixedValue':
                j = NameBCs.index(i)
                for k in range(int(startFace[j]), int(startFace[j])+int(nFaces[j])):
                    o = owner[k]
                    OutputFaces[k,dim] = valor[dim]
                    DeslocamentoBC[k,dim] = valor[dim]
    
    for i in range(Npoints): #interpolando os resultados
        k=0
        for j in range(Nfaces):
            if i in faces[j]:
                DeslocamentoA[i,dim]+=OutputFaces[j,dim]
                k+=1
        DeslocamentoV[i,dim]=DeslocamentoA[i,dim]/k
        
    for i in range(Nfaces): #corrigindo os resultados das BC
        if not np.isnan(DeslocamentoBC[i,dim]):
            for j in range(len(faces[0])):        
                DeslocamentoV[int(faces[i][j]),dim] = DeslocamentoBC[i,dim]
            
    return (DeslocamentoV[:,dim])

#=============================================================================#

def quality():
    dVet = np.zeros((Nneighbour,3))
    ortho = np.zeros(Nneighbour); avg_ortho = np.zeros(Nvolumes); min_ortho = avg_ortho
    owner2 = np.zeros(Nneighbour); qtdFacesVol = np.zeros(Nvolumes)
    for i in range(Nneighbour): # loop faces internas
        o = owner[i] ; n = neighbour[i]
        owner2[i] = owner[i]
        dVet[i] = cVol[n] - cVol[o] # distancia entre centros vol P e N
        ortho[i]=90-np.rad2deg(np.arccos(round(np.dot(dVet[i], areaFaceV[i])/(np.linalg.norm(dVet[i])*np.linalg.norm(areaFaceV[i])),12)))
    for i in range(Nvolumes): # loop em todos os volumes
        # identifica as faces que pertencem ao volume (i)
        facesofvolumeO = np.where(owner2==i)[0]
        # idem so que na lista neighbour pois face nomeada uma unica vez nas listas
        facesofvolumeN = np.where(neighbour==i)[0]
        # quantidade de faces em cada volume
        qtdFacesVol[i] = qtdFacesVol[i] + len(facesofvolumeO) + len(facesofvolumeN)
        # ortogonalidade média no volume
        avg_ortho[i] = (np.sum(ortho[facesofvolumeO],0) + np.sum(ortho[facesofvolumeN],0))/qtdFacesVol[i]
        # mínima ortogonalidade no volume
        if (len(facesofvolumeN) == 0):
            min_ortho[i] = np.min(ortho[facesofvolumeO])
        elif (len(facesofvolumeO) == 0):
           min_ortho[i] = np.min(ortho[facesofvolumeN])
        else:
           min_ortho_prog = np.min(ortho[facesofvolumeO]),np.min(ortho[facesofvolumeN])
           min_ortho[i] = np.min(min_ortho_prog)
    return (avg_ortho, min_ortho)

#=============================================================================#
def write_scalar(dirFinal,varName,solution): # func que escreve resultados

    pathFinal = './' + dirFinal + '/' + varName
    if not os.path.exists(dirFinal):
        os.makedirs(dirFinal)
    
    res = str(solution).replace(']', ');')
    #correção de espaços, parentesis e pontos    
    res = res.replace('[','(')
    res = res.replace('((','(')
    res = res.replace('))',')')
    res = res.replace(' (','(')
    res = res.replace(' )',')')
    res = res.replace(' )',')')
    res = res.replace(' )',')')
    res = res.replace('( ','(')
    res = res.replace('( ','(')
    res = res.replace('.)',')')
    res = res.replace('  ',' ')
    res = res.replace('  ',' ')
    res = res.replace('  ',' ')
    res = ''.join(res)
    
    header = '/*--------------------------------*- C++ -*----------------------------------*\\'
    header += '\n| =========                 |                                                 |'
    header += '\n| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'
    header += '\n|  \\\    /   O peration     | Version:  5.0                                   |'
    header += '\n|   \\\  /    A nd           | Web:      www.OpenFOAM.org                      |'
    header += '\n|    \\\/     M anipulation  |                                                 |'
    header += '\n\*---------------------------------------------------------------------------*/'
    header += '\nFoamFile'
    header += '\n{'
    header += '\n    version     2.0;'
    header += '\n    format      ascii;'
    header += '\n    class       volScalarField;'
    header += '\n    location    "%s";' % dirFinal
    header += '\n    object      %s;' % varName
    header += '\n}'
    header += '\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
    header += '\n'
    header += '\ndimensions      [0 0 0 0 0 0 0];'
    header += '\n'
    header += '\ninternalField nonuniform List<scalar>\n%d' % Nvolumes
    #header += '\n\n'
    
    end = '\n'
    end += '\nboundaryField'
    end += '\n{'
    for i in range(len(NameBCs)):
        end += '\n\t%s' % str(NameBCs[i])
        end += '\n\t{'
        end += '\n\t\ttype fixedValue;'
        end += '\n\t\tvalue uniform 0;'    
        end += '\n\t}'
    end += '\n}'
    end += '\n'
    end += '\n'
    end += '\n// ************************************************************************* //'

    with open(pathFinal, 'w') as outfile:
            outfile.write(header)
            for item in res:            
                outfile.write(item)
            outfile.write(end)

    return ()
###############################################################################
# ------------------------------ SOLVER --------------------------------------#
###############################################################################

intro =  '\n================================================================='
intro += '\n================================================================='
intro += '\n----------------------------- MovMesh ---------------------------'
intro += '\n================================================================='
intro += '\n================================================================='
print (intro)

StartTimeTot = timer() # contabilizador do tempo

#==============================================================================
# -------------------------- CORREÇÕES DOS INPUTS ----------------------------#
#==============================================================================
tnow = tinitial + dt
dirInicial = str(tinitial)
dirFinalr = str(tnow) + '/polyMesh'
dirMalha = 'constant/polyMesh'
dirInicial ='./' + dirInicial + '/' + CI
dirFinal = './' + dirFinalr + '/'
dirMalhaBoundary = './' + dirMalha +'/boundary'
dirMalhaPoints = './' + dirMalha + '/points'
dirMalhaFaces = './' + dirMalha +'/faces'
dirMalhaOwner = './' + dirMalha +'/owner'
dirMalhaNeighbour = './' + dirMalha +'/neighbour'
#=============================================================================#
while tinitial < tfinal:
    ### LEITURA DA MALHA E DAS CONDIÇÕES DE CONTORNO
    ## -- Computa nós da malha
    print('TIMESTEP = %f' % tnow)
    print('Lendo arquivo de malha...')
    timeMeshB = timer()
    points, Npoints = read_file(dirMalhaPoints)
    points = np.asarray(points)
    points = np.array(points, dtype=float)
    #=============================================================================#
    ## -- Computa faces da malha
    faces, Nfaces = read_file(dirMalhaFaces)
    
    for i in range(len(faces)): # necessario as vezes...
        faces[i] = np.array(faces[i])
        faces[i] = faces[i].astype(np.int)
    #=============================================================================#
    ## -- Computa volumes donos das faces na malha
    owner = read_scalarlist(dirMalhaOwner)
    owner = np.asarray(owner);owner = owner.astype(np.int)
    Nowner = len(owner)
    #=============================================================================#
    ## -- Computa volumes vizinhos das faces na malha
    neighbour = read_scalarlist(dirMalhaNeighbour)
    neighbour = np.asarray(neighbour);neighbour = neighbour.astype(np.int)
    Nneighbour = len(neighbour)
    #=============================================================================#
    ## -- Computa a area (modulo e vetor) de cada face e os centros
    areaFace, areaFaceV, cFace = calc_mesh_faces()
    #tempo gasto
    #=============================================================================#
    ## -- Computa o numero de volumes da malha e as coords centroides dos vols
    Nvolumes, cVol, Vol = calc_mesh_vol()
    #tempo gasto
    timeMeshD = timer()
    WallTimeMesh = timeMeshD - timeMeshB
    print ('Leitura da malha terminada em %f segundos.' % WallTimeMesh)
    print ('-----------------------------------------------------------------')
    ## -- Leitura das condições de contorno
    print('Lendo condições de contorno...')
    timeBCB = timer()
    ## -- Computa infos de condicoes de contorno
    NameBCs, TipoBCs, nFaces, startFace = read_bcs_mesh()
    timeBCD = timer()
    BCTime = timeBCD - timeBCB
    print ('Leitura das condições de contorno terminada em %f segundos.' % BCTime)
    print ('-----------------------------------------------------------------')
    #=============================================================================#
    ## -- Cálculo da qualidade da malha
    print ('Calculando métricas de qualidade da malha inicial...')
    QualityB = timer()
    dirFinalscalar = str(tinitial)
    avg_ortho, min_ortho = quality()
    write_scalar(dirFinalscalar, scalar, min_ortho)
    QualityD = timer()
    QualityTime = QualityD - QualityB
    print ('Métricas calculadas em %f segundos.' % QualityTime)
    print ('-----------------------------------------------------------------')
    #=============================================================================#
    #=============================================================================#
    ### SOLUÇÃO DA EQUAÇÃO DE DIFUSÃO
    print('Resolvendo o sistema linear...')
    timeSolB = timer()
    ResultadoFinal = np.zeros([Npoints,3])#inicializando o vetor resultado
    StartTimeSisLin = timer()
    for dim in range(3): #loop em 3D
        A, Su = solucao() #matriz dos coeficientes
        Solucao = bicg(A, Su) #Solucao = bicgstab(A, Su) #alternativa
        Output = Solucao[0] #ignorando uma parte da solução
        ResultadoFinal[:,dim]=interpolator_corrector() #interpolando e corrigindo o resultado
    timeSolD = timer()
    timeSol = timeSolD - timeSolB
    print ('Sistema linear resolvido em %f segundos.' % timeSol)
    print ('-----------------------------------------------------------------')
    #=============================================================================#
    ## -- Escrita do arquivo de saída
    print ('Exportanto resultados...')
    ExportB = timer()
    exportar = points + ResultadoFinal
    write_OF(dirFinal,'points', exportar)      
    ExportD = timer()
    ExportTime = ExportD - ExportB
    print ('Resultados exportados em %f segundos.' % ExportTime)
    print ('-----------------------------------------------------------------')
    #=============================================================================#
    #=============================================================================#
    ### ATUALIZAÇÃO DAS REFERÊNCIAS
    tinitial = tnow    
    tnow = tinitial + dt
    dirInicial = str(tinitial)
    dirFinalr = str(tnow) + '/polyMesh'
    dirMalha = 'constant/polyMesh'
    dirInicial ='./' + dirInicial + '/' + CI
    dirFinal = './' + dirFinalr + '/'
    dirMalhaBoundary = './' + dirMalha +'/boundary'
    dirMalhaPoints = './' + str(tinitial) + '/polyMesh/points'
    dirMalhaFaces = './' + dirMalha +'/faces'
    dirMalhaOwner = './' + dirMalha +'/owner'
    dirMalhaNeighbour = './' + dirMalha +'/neighbour'
    #=============================================================================#
    #=============================================================================#
    #fim do while
### LEITURA DA MALHA E CÁLCULOS DA MALHA
## -- Computa nós da malha
print('Lendo último arquivo de malha...')
timeMeshB = timer()
points, Npoints = read_file(dirMalhaPoints)
points = np.asarray(points);points = points.astype(np.float)
#=============================================================================#
## -- Computa faces da malha
faces, Nfaces = read_file(dirMalhaFaces)  
for i in range(len(faces)): # necessario as vezes...
    faces[i] = np.array(faces[i])
    faces[i] = faces[i].astype(np.int)
#=============================================================================#
## -- Computa volumes donos das faces na malha
owner = read_scalarlist(dirMalhaOwner)
owner = np.asarray(owner);owner = owner.astype(np.int)
Nowner = len(owner)
#=============================================================================#
## -- Computa volumes vizinhos das faces na malha
neighbour = read_scalarlist(dirMalhaNeighbour)
neighbour = np.asarray(neighbour);neighbour = neighbour.astype(np.int)
Nneighbour = len(neighbour)
#=============================================================================#
## -- Computa a area (modulo e vetor) de cada face e os centros
areaFace, areaFaceV, cFace = calc_mesh_faces()
#tempo gasto
#=============================================================================#
## -- Computa o numero de volumes da malha e as coords centroides dos vols
Nvolumes, cVol, Vol = calc_mesh_vol()
#tempo gasto
timeMeshD = timer()
WallTimeMesh = timeMeshD - timeMeshB
##=============================================================================#
## -- Cálculo da qualidade da malha
print ('Calculando métricas de qualidade da malha...')
QualityB = timer()
dirFinalscalar = str(tinitial)
avg_ortho, min_ortho = quality()
write_scalar(dirFinalscalar, scalar, min_ortho)
QualityD = timer()
QualityTime = QualityD - QualityB
print ('Métricas calculadas em %f segundos.' % QualityTime)
print ('-----------------------------------------------------------------')
#=============================================================================#
#=============================================================================#
### FIM DO CÓDIGO
TimeFinal = timer()
WallTimeTot = TimeFinal - StartTimeTot
print ('Código finalizado em %f segundos.' % WallTimeTot)
fim =  '================================================================='
fim += '\n================================================================='
fim += '\n----------------------- Código finalizado -----------------------'
fim += '\n================================================================='
fim += '\n================================================================='
print (fim)
#fim do código
#=============================================================================#
#=============================================================================#
#=============================================================================#
#=============================================================================#
