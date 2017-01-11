import ctypes
from ctypes.util import find_library
from numpy.ctypeslib import ndpointer



# find and load the library
lib_spec = ctypes.cdll.LoadLibrary(find_library('spec'))


# set the argument type for diagonalization oneliner
lib_spec.spec_pencil.argtypes = [ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), 
								 ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), 
								 ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), 
								 ctypes.c_int]

lib_spec.spec_pencil.restype = None


def spec_pencil(eig_val, eig_vect, hreal, himag, sreal, simag, n):
	return lib_spec.spec_pencil(eig_val, eig_vect, hreal, himag, sreal, simag, n)
