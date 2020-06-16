#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

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


## -- Main

#respiracao T=3.75s
dt = 0.05
# - Respiracao 1
tr1= np.arange(0, 15.05, dt)
r1 = -3*np.sin((tr1/3.75)*2*np.pi)
resp1 = np.column_stack((tr1, r1))

# - Respiracao 2
tr2= np.arange(15.4, 30.45, dt)
r2 = -3*np.sin((tr1/3.75)*2*np.pi)
resp2 = np.column_stack((tr2, r2))

# Tosse
tt = np.arange(0, 0.41, 0.01)

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

tt = tt + 15
tosse = np.column_stack((tt, tos))


# - Plots pra verificacao
plt.figure(figsize=(6.4,5),dpi=200)

plt.plot(tr1,r1,'r',label='Resp 1')
plt.plot(tr2,r2,'b',label='Resp 2')
plt.plot(tt,tos,'k',label='Tosse')

plt.xlabel('Time [t(s)]')
plt.ylabel('r(1/2)')
plt.legend()
plt.tight_layout(pad=0.5)


plt.figure(figsize=(6.4,5),dpi=200)

plt.plot(tt,tos,'k',label='Tosse')

plt.xlabel('Time [t(s)]')
plt.ylabel('r(1/2)')
plt.legend()
plt.tight_layout(pad=0.5)

plt.show()


# - Escreve arquivos entrada converge
HEADER = 'TEMPORAL\nSEQUENTIAL\nsecond\tAverage_velocity'
np.savetxt('resp1.in',resp1,header=HEADER,delimiter='\t',fmt=['%.2f','%.9f'],comments='')
np.savetxt('resp2.in',resp2,header=HEADER,delimiter='\t',fmt=['%.2f','%.9f'],comments='')
np.savetxt('tosse.in',tosse,header=HEADER,delimiter='\t',fmt=['%.2f','%.9f'],comments='')