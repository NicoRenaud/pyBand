import numpy as np
from string import rjust
from string import ljust
import commands
import os
import sys
import copy

################################################
#
#	print a mopac file for MO
#
################################################
def print_mo_mopac(wH,uH,s,indexHomo,system,filename,PARAM,nb_print_orb,kpt):

	natom = system.nAtom
	natom_print = natom

	#nb_orb = len(wH)
	nb_orb = len(uH)
	count = 0

	# comute the inverse of S
	invS = np.linalg.inv(s)


	# open the output file
	f = open(filename,'w')

	# header
	f.write("        %d MOPAC-Graphical data Version 2012.13.084W K = [%1.3f %1.3f %1.3f]\n" %(natom_print,kpt[0],kpt[1],kpt[2]))

	# print the atoms
	for i in range(natom):
		at = find_index(system.atomType[i])
		f.write('%4d %*f%*f%*f  0.0000\n' %(at,12,system.xyz[i][0],12,system.xyz[i][1],12,system.xyz[i][2]))


	# print the slater exponents
	for i in range(natom):
		sc1,sc2,sc3 = get_slater_coeff(system.atomType[i],PARAM)
		f.write("  %1.7f  %1.7f  %1.7f\n" %(sc1,sc2,sc3))

	# determine what to print
	index_print_orb = []
	occ = []

	for i in range(nb_print_orb):
		index_print_orb.append(indexHomo-nb_print_orb/2+1+i)
		if i<nb_print_orb/2:
			occ.append(2)
		else:
			occ.append(0)

	# print the orbitals
	for iprint in range(nb_print_orb):

		iorb = index_print_orb[iprint]
		f.write(" ORBITAL %3d  A  %2.5g\n" %(occ[iprint],wH[iorb]-wH[indexHomo]))

		for jcomp in range(nb_orb):
			f.write("% 1.8E" %(uH[jcomp,iorb]).real)
			count += 1
			if count == 5:
				f.write("\n")
				count  = 0
		if count>0:
				f.write("\n")
				count  = 0

	# print the inverse matrix
	count = 0
	f.write("INVERSE_MATRIX[%dx%d]=\n"%(nb_orb,nb_orb))
	for i in range(nb_orb):
		for j in range(i+1):
			f.write("% 1.8E" %(invS[j,i]))
			count+=1
			if count == 5:
				f.write("\n")
				count  = 0
		if count>0:
				f.write("\n")
				count  = 0
	f.write(" Keywords: SYMMETRY GRAPHF")
	f.close()


###########################################################
#
# print the DOS
#
###########################################################
def print_dos(w,index_homo,filename,width,emin,emax,nE):

	f = open(filename,'w')
	E = np.linspace(emin,emax,nE)
	for iE in range(len(E)):
		g = 0
		e = E[iE]
		for iW in range(len(w)):
			w0 = w[iW]-w[index_homo]
			g += np.exp( -(e-w0)**2/width**2 )
		f.write("%f %f\n" %(e,g))
	f.close()

###########################################################
#
# find the index that corresponds to a given atom name
#
###########################################################
def find_index(name):
	#list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
	# replaced Fe by Pb otherwise Jmol can't print MO on Pb ....
	list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "PB", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Fe", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
	index= [i for i, x in enumerate(list) if x == name]
	return index[0]+1



###########################################################
#
# find the slater expoants of a given atom
#
###########################################################
def get_slater_coeff(name,param):

	f =open(param,'r')
	data = f.readlines()
	f.close

	for i in range(len(data)):
		data[i] = data[i].split()
		if data[i][0] == name:
			sc1 = data[i][7]
			sc2 = data[i][13]
			sc3 = data[i][19]
			break
	return float(sc1),float(sc2),float(sc3)



###########################################################
#
# Get the number of atomic orbitals per atoms
#
###########################################################
def get_nbAO(atomType,param):

	nbAo = []

	f =open(param,'r')
	data = f.readlines()
	f.close

	for i in range(len(data)):
		data[i] = data[i].split()

	for i in range(len(atomType)):

		for j in range(len(data)):
			if data[j][0] == atomType[i]:
				if int(data[j][2]) != 0:
					nbAo.append([i])
				if int(data[j][3]) != 0:
					nbAo.append([i,i,i])
				if int(data[j][4]) != 0:
					print 'd orbitals detected '
					nbAo.append([i,i,i,i,i])
				break;
	nbAo = sum(nbAo,[])
	return nbAo
