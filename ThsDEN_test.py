from estimator.lwe_parameters import LWEParameters as LWEest
from estimator import *
import math
etas,etae=2,2
Q=7681
# Q=10753
n=256
d=3 # d=5 for 256-bit security
N=n*d
params = LWEest(
     n=N,
     q=Q,#1280,256 bit
     Xs=ND.DiscreteGaussian(math.sqrt(etas/2)), # s
     Xe=ND.DiscreteGaussian(math.sqrt(etae/2)), # e
     m=N,
     tag="example"
)

r=LWE.estimate(params)
