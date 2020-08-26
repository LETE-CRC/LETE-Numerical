#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
                 Python code for elevator time LETE-(no)COVID
   Created by Combustion Research Center CRC at LETE - Sao Paulo, Brasil
   Laboratory of Environmental and Thermal Engineering - LETE
   Escola Politecnica da USP - EPUSP

===============================================================================
version:0.0 - 06/2020: Helio Villanueva
"""

import numpy as np
import matplotlib.pyplot as plt

Vel_Nom = 1 # [m/s] p/ elevadores de 6 a 8 pessoas
Accel_time = 3 # [s]
Door_time = 3.9 # [s]
Travel_time = 10 #[s]
t_resp = 3.75 # [s] duracao respiracao
t_tosse = 0.41 # [s] duracao tosse
t_inicio_tosse = 5.7 # [s]
dt = 0.05
pos_fechada = 0.48 #[m] dist percorrida p fechar uma folha da porta

# -----------------------------------------------------------------------------
# Definicao do tempo
tf = 2*Accel_time + 2*Door_time + Travel_time
time = np.arange(0,tf+dt,dt)

# -----------------------------------------------------------------------------
# - aceleracao do elevador
a = Vel_Nom/Accel_time
accel = np.zeros_like(time)
accel[(time > Door_time) & (time < Door_time+Accel_time)] = a
accel[(time > Travel_time+Door_time+Accel_time) & (time < time[-1]-Door_time)] = -a

# -----------------------------------------------------------------------------
# - velocidade do elevador
vel = np.zeros_like(time)
for i,v in enumerate(vel):
    if i==0:
        vel[i]=0
    else:
        vel[i] = vel[i-1]+accel[i]*dt

# -----------------------------------------------------------------------------
# - artificial g
''' Supondo gravidade no sentido -Y
    Caso elevador subindo: - accel
    Caso elevador descendo: + accel '''
gl = -9.81 + accel
gravity = np.column_stack((time, gl))

# -----------------------------------------------------------------------------
# - Respiracao 1
tr1= np.arange(0, t_inicio_tosse, dt)
r1 = -3*np.sin((tr1/t_resp)*2*np.pi)
r1[-1] = 0
resp1 = np.column_stack((tr1, r1))

# - Respiracao 2
#tr2= np.arange(t_inicio_tosse+t_tosse, tf, dt)
tr2= np.arange(0, tf-(t_inicio_tosse+t_tosse), dt)
r2 = -3*np.sin((tr2/t_resp)*2*np.pi)
tr2 = tr2 + t_inicio_tosse + t_tosse
resp2 = np.column_stack((tr2, r2))

Resp = np.append(resp1,resp2,axis=0)

# - Tosse
tt = np.arange(0, t_tosse, 0.01)

coeffs = [18693663.2299235,
          -32465014.3862506,
          23226581.353063,
          -8783010.36337336,
          1861672.52173274,
          -212404.400366862,
          10580.039007484,
          -40.1864550037,
          0.4068601721]

tos = np.polyval(coeffs,tt)
tos[0] = 0
tos[-1] = 0

tt = tt + t_inicio_tosse
tosse = np.column_stack((tt, tos))
tosse = np.append([[0,0]],tosse,axis=0)
tosse = np.append(tosse,[[tf,0]],axis=0)

# -----------------------------------------------------------------------------
# - movimentacao das portas
p_pos = np.zeros_like(time)
for i,t in enumerate(time):
    if i==0:
        p_pos[i] = 0
    elif t<Door_time:
        p_pos[i] = p_pos[i-1] + dt*pos_fechada/Door_time
    elif t>Door_time+2*Accel_time+Travel_time:
        p_pos[i] = p_pos[i-1] - dt*pos_fechada/Door_time
    else:
        p_pos[i] = pos_fechada

portas = np.column_stack((time, p_pos))

# -----------------------------------------------------------------------------
# - salva tabela no formato do converge
HEADERg = 'TEMPORAL\nSEQUENTIAL\nsecond\tgravity_y'
HEADERpos = 'TEMPORAL\nSEQUENTIAL\nsecond\tlift'
HEADER = 'TEMPORAL\nSEQUENTIAL\nsecond\tAverage_velocity'
np.savetxt('gravity.in',gravity,header=HEADERg,delimiter='\t',fmt=['%.2f','%.9f'],comments='')
np.savetxt('resp1.in',resp1,header=HEADER,delimiter='\t',fmt=['%.2f','%.9f'],comments='')
np.savetxt('resp2.in',resp2,header=HEADER,delimiter='\t',fmt=['%.2f','%.9f'],comments='')
np.savetxt('resp.in',Resp,header=HEADER,delimiter='\t',fmt=['%.2f','%.9f'],comments='')
np.savetxt('tosse.in',tosse,header=HEADER,delimiter='\t',fmt=['%.2f','%.9f'],comments='')
np.savetxt('portas.in',portas,header=HEADERpos,delimiter='\t',fmt=['%.2f','%.9f'],comments='')

# -----------------------------------------------------------------------------
# - Def visual plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes',linewidth=2,labelsize=18)
plt.rc('axes.spines',top=0,right=0)

plt.rc('xtick',labelsize=16)
plt.rc('xtick.major',size=5,width=2)
plt.rc('xtick.minor',visible=1)
plt.rc('ytick',labelsize=16)
plt.rc('ytick.major',size=5,width=2)
plt.rc('ytick.minor',visible=1)

# - Plots
plt.figure(figsize=(10,5),dpi=200)
ax = plt.subplot(321)
plt.plot(time,vel,'k')
plt.ylabel('V $[m/s]$')
plt.ylim(0,1.01)
plt.grid(which='minor',linestyle='--')
plt.grid(which='major',linestyle='--')

plt.subplot(323,sharex=ax)
plt.plot(time,accel,'k')
plt.ylabel('a $[m/s^2]$')
plt.grid(which='minor',linestyle='--')
plt.grid(which='major',linestyle='--')

plt.subplot(325,sharex=ax)
plt.plot(time,gl,'k')
plt.ylabel("g' $[m/s^2]$")
plt.grid(which='minor',linestyle='--')
plt.grid(which='major',linestyle='--')

#plt.xlabel('tempo [s]')
#plt.xlim(time[0],time[-1])
#plt.tight_layout(pad=0.5)

#plt.figure(figsize=(6.4,5),dpi=200)
plt.subplot(322,sharex=ax)
plt.plot(tr1,r1,'r',label='Resp 1')
plt.plot(tr2,r2,'b',label='Resp 2')
plt.plot(tt,tos,'k',label='Tosse')
plt.ylabel('v $[m/s]$')
plt.legend()

plt.subplot(324,sharex=ax)
plt.plot(time,p_pos,'k')
plt.ylabel('portas $[m]$')
plt.ylim(0,0.501)
plt.grid(which='minor',linestyle='--')
plt.grid(which='major',linestyle='--')

plt.xlabel('tempo [s]')
plt.xlim(time[0],time[-1])
plt.tight_layout(pad=0.5)
plt.show()