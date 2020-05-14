import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
#respiracao T=3.75s
tr1= np.arange(0, 15.05, 0.05).tolist()
tr1 = np.transpose(tr1)
r1 = -3*np.sin((tr1/3.75)*2*np.pi)
resp1 = np.column_stack((tr1, r1))
tr2= np.arange(15.4, 30.45, 0.05).tolist()
tr2 = np.transpose(tr2)
r2 = -3*np.sin((tr1/3.75)*2*np.pi)
resp2 = np.column_stack((tr2, r2))
#tosse
tt = np.arange(0, 0.41, 0.01).tolist()
tt = np.transpose(tt)
c1=18693663.2299235
c2=-32465014.3862506
c3=23226581.353063
c4=-8783010.36337336
c5=1861672.52173274
c6=-212404.400366862
c7=10580.039007484
c8=-40.1864550037
c9=0.4068601721
tos = c1*tt**8 + c2*tt**7 + c3*tt**6 + c4*tt**5 + c5*tt**4 + c6*tt**3 + c7*tt**2 + c8*tt + c9
tt = tt + 15
tosse = np.column_stack((tt, tos))

print(resp2)
