#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


#respiracao T=3.75s
dt = 0.05
# - Respiracao 1
tr1= np.arange(0, 15.05, dt)
r1 = -3*np.sin((tr1/3.75)*2*np.pi)
resp1 = np.column_stack((tr1, r1))

# - Respiracao 1
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

print(resp2)
