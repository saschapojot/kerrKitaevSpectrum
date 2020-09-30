import numpy as np
#fermion number
N=500
dk=2*np.pi/N
kIndHalf=range(0,int(N))

#parameters before quench
mu0=0
t0=1.0
d0=-1.0

#parameters after the quench
mu1=1
t1=t0
d1=d0
lmd=3
lmdAll=range(0,20)
threadNum=12
#occupation ratio
rho=1

tol=1e-15

#dict of spectrum
deltaKAll=dict()
EkAll=dict()
linEAll=[]

h0Val = lmd / 2 * rho - mu1 / 2

# solve deltaK, coefs of 4th order algebraic equation

coef4Val = lmd ** 2 / 4
#linear spectrum
h0LVal=-mu1/2
