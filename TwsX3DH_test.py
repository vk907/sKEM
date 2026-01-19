from estimator.lwe_parameters import LWEParameters as LWEest
from estimator import *
import math
etas,etae=3,4
# etas,etae=2,2 for 256-bit security
Q=7681
n=256
d=2
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
