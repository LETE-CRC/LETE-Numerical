#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 2015 by helio
Modified on Tue Oct 24 2017 by filipi
"""

import numpy as np
import re
from scipy.sparse.linalg import bicg # bicgstab
import shutil
import os
from timeit import default_timer as timer

#==============================================================================
# --------------------- INPUTS E SELEÇÕES NECESSÁRIAS ------------------------#
#==============================================================================

## -- Propriedades
gamma = 0.1 # RIGIDEZ DA MALHA

LaplacianScheme = 'linear'


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
    with open('./constant/polyMesh/boundary', 'r') as infile: # abre o arquivo a ser lido
        arquivo = infile.read() # le o arquivo como unica string
        arquivo = ''.join(arquivo.split()) # retira todos espacos em branco
        namebcs = re.findall(r'(\w+)\{type.\}*', arquivo) # encontra nomes das bcs
        tipobc = re.findall(r'type(\w+)', arquivo) # encontra tipo
        nfaces = re.findall(r'nFaces(\w+)', arquivo) # encontra id primeira face na lista de faces
        startface = re.findall(r'startFace(\w+)', arquivo) # encontra id primeira face na lista de faces
        nfaces = np.asarray(nfaces);nfaces = nfaces.astype(np.int)
    return (namebcs, tipobc, nfaces, startface)

#=============================================================================#

def write_OF(varName,solution): # func que escreve resultados
    '''
    Write OpenFOAM format solution
        varName -> the name of the variable to write
        ex: 'T'
        
        solution -> the variable
        ex: T
    '''
    np.set_printoptions(threshold='nan')
    dirInicial = './0/'
    dirFinal = './1/'
    pathInicial = dirInicial + varName
    pathFinal = dirFinal + varName
    if not os.path.exists(dirFinal):
        os.makedirs(dirFinal)
        
    shutil.copyfile(pathInicial, pathFinal)
    with open(pathInicial, 'r') as f:
        data = f.readlines()

    res = str(solution).replace(']', ');')
    res = res.replace('[','(')
    data.insert(20, res)
    data = "".join(data)
    with open(pathFinal, 'w') as f:
        f.write(data)

    chop = re.compile('internalField\s+uniform\s+0\;')
    with open(pathFinal, 'r') as f:
        data = f.read()
    
    # chop text between #chop-begin and #chop-end
    data_chopped = chop.sub('internalField nonuniform List<scalar>', data)
    # save result
    with open(pathFinal, 'w') as f:
        f.write(data_chopped)

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
    
    DeslocamentoBC[:]=np.nan

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

    with open('./0/T', 'r') as infile:
        arquivo = infile.read()
        arquivo = ''.join(arquivo.split()) # retira todos espacos em branco

    for i in NameBCs:
        my_regex1 = re.escape(i) + r"\{type(\w+)\}*"
        my_regex2 = re.escape(i) + r"\{type\w+\;valueuniform(\w+)\}*"
        TipobcDic = re.findall(my_regex1, arquivo)
        TipobcDic = str(TipobcDic).strip('[\']')
        valor = re.findall(my_regex2, arquivo, re.MULTILINE)
        valor = str(valor).strip('[\']')
        try:
            valor = float(valor)
        except:
            pass
            
        print re.escape(i),'--->', TipobcDic, '--->',valor # Mostra no terminal as condicoes de contorno
        
        if TipobcDic == 'empty' or TipobcDic == 'symmetryPlane':
            j = NameBCs.index(i)
            for k in range(int(startFace[j]), int(startFace[j])+int(nFaces[j])):
                SuBC[k] = 0
                SpBC[k] = 0
            
        if TipobcDic == 'fixedValue':
            j = NameBCs.index(i)
            for k in range(int(startFace[j]), int(startFace[j])+int(nFaces[j])):
                SuBC[k] = (areaFace[k] * gamma / dPf[k]) * valor
                SpBC[k] = -areaFace[k] * gamma / dPf[k]

        
    for i in range(Nfaces):
        o = owner[i] # para cada face, o é o valor da lista owner (n do vol de controle)
        # Su (lado direito) e o inicial (com as geracoes volumetricas) mais contribuicoes cada face
        Su[o] = Su[o] + SuBC[i]
        Sp[o] = Sp[o] + SpBC[i]

    for i in range(len(A)):
        A[i,i] = -np.sum(A[i])-Sp[i] # aP da soma dos vizinhos
        
    return (A, Su, OutputFacesC, DeslocamentoBC)

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
            OutputFaces[i,0] = (Output[o]*dVolO[i]+Output[n]*dVolN[i])/dCent[i] #interpolando o resultado para as faces
    
    with open('./0/T', 'r') as infile:
        arquivo = infile.read()
        arquivo = ''.join(arquivo.split()) # retira todos espacos em branco
    
    for i in NameBCs:
            my_regex1 = re.escape(i) + r"\{type(\w+)\}*"
            my_regex2 = re.escape(i) + r"\{type\w+\;valueuniform(\w+)\}*"
            TipobcDic = re.findall(my_regex1, arquivo)
            TipobcDic = str(TipobcDic).strip('[\']')
            valor = re.findall(my_regex2, arquivo, re.MULTILINE)
            valor = str(valor).strip('[\']')
            try:
                valor = float(valor)
            except:
                pass
                
            print re.escape(i),'--->', TipobcDic, '--->',valor # Mostra no terminal as condicoes de contorno
            
            if TipobcDic == 'empty' or TipobcDic == 'symmetryPlane':
                j = NameBCs.index(i)
                for k in range(int(startFace[j]), int(startFace[j])+int(nFaces[j])):
                    o = owner[k]
                    OutputFaces[k,0] = Output[o]             
                    
            if TipobcDic == 'fixedValue':
                j = NameBCs.index(i)
                for k in range(int(startFace[j]), int(startFace[j])+int(nFaces[j])):
                    o = owner[k]
                    OutputFaces[k,0] = valor
                    DeslocamentoBC[k,0] = valor
    
    for i in range(Npoints): #interpolando os resultados
        k=0
        for j in range(Nfaces):
            if i in faces[j]:
                DeslocamentoA[i,0]+=OutputFaces[j,0]
                k+=1
        DeslocamentoV[i,0]=DeslocamentoA[i,0]/k
        
    for i in range(Nfaces): #corrigindo os resultados das BC
        if not np.isnan(DeslocamentoBC[i,0]):
            DeslocamentoV[int(faces[i][0]),0] = DeslocamentoBC[i,0]
            DeslocamentoV[int(faces[i][1]),0] = DeslocamentoBC[i,0]
            DeslocamentoV[int(faces[i][2]),0] = DeslocamentoBC[i,0]
            DeslocamentoV[int(faces[i][3]),0] = DeslocamentoBC[i,0]
        
    return (DeslocamentoV)

###############################################################################
# ------------------------------ SOLVER --------------------------------------#
###############################################################################

intro = '\n============================================='
intro += '\n============================================='
intro += '\n------------------ MovMesh ------------------'
intro += '\n============================================='
intro += '\n============================================='
print (intro)

StartTimeTot = timer() # contabilizador do tempo

### LEITURA DA MALHA
## -- Computa nós da malha
points, Npoints = read_file('./constant/polyMesh/points')
points = np.asarray(points);points = points.astype(np.float)
#=============================================================================#
## -- Computa faces da malha
faces, Nfaces = read_file('./constant/polyMesh/faces')

for i in range(len(faces)): # necessario as vezes...
    faces[i] = np.array(faces[i])
    faces[i] = faces[i].astype(np.int)
#=============================================================================#
## -- Computa volumes donos das faces na malha
owner = read_scalarlist('./constant/polyMesh/owner')
owner = np.asarray(owner);owner = owner.astype(np.int)
Nowner = len(owner)
#=============================================================================#
## -- Computa volumes vizinhos das faces na malha
neighbour = read_scalarlist('./constant/polyMesh/neighbour')
neighbour = np.asarray(neighbour);neighbour = neighbour.astype(np.int)
Nneighbour = len(neighbour)
#=============================================================================#
timeFiles = timer()
#=============================================================================#
## -- Computa a area (modulo e vetor) de cada face e os centros
print(' Calculando áreas e centros das faces da malha')
areaFace, areaFaceV, cFace = calc_mesh_faces()
#tempo gasto
timeA = timer()
timeArea = timer() - timeFiles
print("done: %f s " %timeArea)
print(' Calculando volumes e centros dos volumes da malha')
#=============================================================================#
## -- Computa o numero de volumes da malha e as coords centroides dos vols
Nvolumes, cVol, Vol = calc_mesh_vol()
#tempo gasto
timeVol = timer() - timeA
print("done: %f s " %timeVol)
timeMesh = timer()
#=============================================================================#
#=============================================================================#
### CONDIÇÕES DE CONTORNO
print ('----------- Condições de contorno -----------\n')
## -- Computa infos de condicoes de contorno
NameBCs, TipoBCs, nFaces, startFace = read_bcs_mesh()
#=============================================================================#
## -- Computa a matriz dos coeficientes e vetor lado direito
A, Su, OutputFacesC, DeslocamentoBC = solucao()
#=============================================================================#
#=============================================================================#
### SOLUÇÃO DA EQUAÇÃO DE DIFUSÃO
## -- Solução do sistema linear
StartTimeSisLin = timer()
Solucao = bicg(A, Su)
#Solucao = bicgstab(A, Su) #alternativa
Output = Solucao[0]
#tempo gasto
TimeFinal = timer()
WallTimeSisLin = TimeFinal - StartTimeSisLin
WallTimeTot = TimeFinal - StartTimeTot
WallTimeFiles = timeFiles - StartTimeTot
WallTimeMesh = timeMesh - timeFiles
WallTimeAssembly = StartTimeSisLin - timeMesh
#=============================================================================#
## -- Interpolação e correção do resultado
ResultadoFinal=interpolator_corrector()
#=============================================================================#
## -- Impressão no terminal dos resultados para conferência
print ('Volumes: %d' %Nvolumes)
final1 = '\n============================================='
final1 += '\n----------------- Solução  ------------------\n\n'
final2 = '\n============================================='
print (final1 + str(ResultadoFinal) + final2)
#=============================================================================#
## - Impressão do tempo total gasto para a solução final
print ('Tempo para leitura dos arquivos: %f segundos' % WallTimeFiles)
print ('Tempo para cálculos da malha: %f segundos' % WallTimeMesh)
print ('Tempo para cálculos das matrizes: %f segundos' % WallTimeAssembly)
print ('Tempo para solucao sistema linear: %f segundos' % WallTimeSisLin)
print ('Tempo total: %f segundos' % WallTimeTot)
#=============================================================================#
## -- Escrita do arquivo de saída
#write_OF('T',T)
#=============================================================================#
## -- Mensagem de aviso de fim de código
fim = '\n============================================='
fim += '\n============================================='
fim += '\n------------- Código finalizado -------------'
fim += '\n============================================='
fim += '\n============================================='
print (fim)

#fim do código
#=============================================================================#
#=============================================================================#
#=============================================================================#
#=============================================================================#
