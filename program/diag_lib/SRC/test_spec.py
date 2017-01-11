import spec_module as spec
import numpy as np
import ctypes
import scipy.linalg as scila
import sys

# show what's inside the module
dir(spec)

# dim
n=2

h = np.zeros((2,2),dtype='complex')
h[0,1] = 1
h[1,0] = 1

# create a S 
s = np.eye(n,dtype='complex')
s[0,1] = 0.1
s[1,0] = 0.1



alpha, beta, vl, vr, work, info = scila.lapack.zggev(h,s)
print alpha/beta
print vl

w,u = scila.eig(h,s)
print w
print u
sys.exit()

# crete void containers
eig_val = np.zeros(n)
eig_vect = np.zeros((n,n))

# diagonalize
spec.spec_pencil(eig_vect, eig_val, h.real, h.imag, s.real, s.imag, n)
eig_vect = eig_vect.T

print "Eigen value with C"
print eig_val

print "Eigen vectors with C"
print eig_vect


