#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from glob import glob
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

files = glob('*.dat')
files.sort()

plt.figure(figsize=(6.4,5),dpi=200)

for file in files:
    name = file.strip('.dat')
    name = name.strip('0')
    name = name.replace('u','$\mu$')
    with open(file) as f:
        data = np.genfromtxt(f,skip_header=1,delimiter=',')
        plt.plot(data[:,0],data[:,1],'+',label=name)

plt.xlabel('Time after injection [t(s)]')
plt.ylabel('Mean vertical position [$\overline{y}$(m)]')
plt.ylim(0,2.6)
plt.xlim(0.1,100)
plt.xscale('log')
plt.legend()
plt.tight_layout(pad=0.5)
plt.savefig('meanVerticalPos.pdf')
plt.show()