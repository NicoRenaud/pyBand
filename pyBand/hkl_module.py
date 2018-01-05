import ctypes
from ctypes.util import find_library
from numpy.ctypeslib import ndpointer
import os

# find and load the library
path = os.path.dirname(os.path.realpath(__file__))
libname = path + '/hkl.so'
if not os.path.isfile(libname):
	raise FileNotFoundError('%s not found' %libname)

# load the library
lib_hkl = ctypes.cdll.LoadLibrary(libname)

# strin pointer type
ctypes_str = ctypes.c_char_p


# set the argument type for number of orbtitals
lib_hkl.compute_nb_orb.argtypes = [ctypes_str]
lib_hkl.compute_nb_orb.restype = ctypes.c_int


# set the argument type for compute_hamiltonian
lib_hkl.compute_huckel_hamiltonian_general.argtypes = [ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ctypes.c_int, ctypes_str]
lib_hkl.compute_huckel_hamiltonian_general.restype = None


# compute the Hamiltonian with huckel
def hklHam(h,s,n,filename):
	return lib_hkl.compute_huckel_hamiltonian_general(h,s,n,filename)

# compute the number of orbital in the system
def nbOrb(filename):
	return lib_hkl.compute_nb_orb(filename)
