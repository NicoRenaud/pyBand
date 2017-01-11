import numpy as np
from string import rjust
from string import ljust
import commands
import os
import sys
import copy
import scipy.linalg as scila

################################################
#
#	 Count the number of electrons
#
################################################
def count_electron(cell,hkl_param):

		param_raw = []
		atom_type = []
		norb_type = []


		#open the input file, read the info, close the file
		inputFile1=open(hkl_param,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()

		#make a matrix of the atom information
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1
		
		for x in range(lenInput):
					param_raw.append(inputFileLines[x].split())
					atom_type.append([param_raw[x][0] , param_raw[x][1]])
					norb_type.append( [ int(param_raw[x][2]),  int(param_raw[x][3]),  int(param_raw[x][4]), int(param_raw[x][5]) ] )
		nelec = 0
		nb_orb = []
		for i in range(cell.nAtom):
			for j in range(len(atom_type)):
				if cell.atomType[i] == atom_type[j][0]:

					# count number of electrons
					nelec += int(atom_type[j][1])

					# count number of orbitals
					no = 0
					if norb_type[j][0]	!= 0:
						no += 1
					if norb_type[j][1]	!= 0:
						no += 3
					if norb_type[j][2] != 0:
						no += 5
					if  norb_type[j][3] != 0:
						no += 7

					nb_orb.append( no )
		return nelec, nb_orb

################################################
#
#	 Block diagonalize  Hcell Hint
#											Hint  Hcell
################################################
def block_diag(hcell,scell,hint,sint):

	
	# diagonalize the cell hamiltonian
	wcell,ucell = np.linalg.eigh(hcell)
	norb = len(hcell)
	
	# index of the two cell in the total matrices
	ind_cell1 = range(norb)
	ind_cell2 = range(norb,2*norb)
	
	# form the oration matrix
	Urot = np.zeros((2*norb,2*norb))
	Urot[np.ix_(ind_cell1,ind_cell1)] = ucell
	Urot[np.ix_(ind_cell2,ind_cell2)] = ucell



	# form the total Hamiltonian
	Htot = np.zeros((2*norb,2*norb))
	Htot[np.ix_(ind_cell1,ind_cell1)] = hcell
	Htot[np.ix_(ind_cell2,ind_cell2)] = hcell
	Htot[np.ix_(ind_cell1,ind_cell2)] = hint
	Htot[np.ix_(ind_cell2,ind_cell1)] = hint.T

	# form the total overlap matrnp.ix
	Stot = np.zeros((2*norb,2*norb))
	Stot[np.ix_(ind_cell1,ind_cell1)] = scell
	Stot[np.ix_(ind_cell2,ind_cell2)] = scell
	Stot[np.ix_(ind_cell1,ind_cell2)] = sint
	Stot[np.ix_(ind_cell2,ind_cell1)] = sint.T

	# rotate the matrices
	Htot = np.dot(Urot.T,np.dot(Htot,Urot))
	Stot = np.dot(Urot.T,np.dot(Stot,Urot))

	# extract the new matrices
	hcell = Htot[np.ix_(ind_cell1,ind_cell1)]
	hint  = Htot[np.ix_(ind_cell1,ind_cell2)]
	scell = Stot[np.ix_(ind_cell1,ind_cell1)]
	sint  = Stot[np.ix_(ind_cell1,ind_cell2)]


	return hcell,scell,hint,sint

################################################
#
#	Orthogonalize the basis
#
################################################


def ortho(H,S):

	n = S.shape[0]
	s, U = eigh(S)
	V = dot(dot(U, diag(1 / sqrt(s))), U.T)
	Ho = dot(dot(V, H), V)

	return Ho

################################################
#
#	Compute the electronic struture
#
################################################

def compute_h(cell1,vect,huckel_exe,hkl_param,keht):

	# name f the inpyt file for huckel
	hkl_fn = "two_cell.in"

	# copyt the cell and translate it
	cell2=copy.deepcopy(cell1)
	cell2.translate(vect)

	# create a input file for huckel
	export_huckel(hkl_fn,cell1,cell2,hkl_param,keht)
	
	# exectute huckel
	cmd = huckel_exe + ' ' + hkl_fn + ' ./'
	os.system(cmd)

	# import the results in pyhton
	hcell = np.loadtxt("hcell.dat")
	scell = np.loadtxt("scell.dat")
	hint = np.loadtxt("hint.dat")
	sint = np.loadtxt("sint.dat")
	
	# save the individual data
	#os.system("rm -f two_cell.in hcell.dat scell.dat hint.dat sint.dat H.dat S.dat")
		
	return hcell, scell, hint, sint


################################################
#
#	Export the husky file for the junction
#
################################################

def export_huckel_old(huckelFileName,cell1,cell2,hkl_param,keht):
	

	natom_cell1 = len(cell1.xyz)
	natom_cell2 = len(cell2.xyz)

	
	
	fp = open(huckelFileName,'w')
	fp.write(" //////////////////////////////\n //   pos for Huckel \n //////////////////////////////\n")
	fp.write("nb_atom %d\n" %(natom_cell1+natom_cell2) )
	fp.write("nb_atom_cell 1 %d\n" %(natom_cell1) )
	fp.write("nb_atom_cell 2 %d\n" %(natom_cell2) )
	fp.write("parameters %s\n" %(hkl_param))
	fp.write("Keht %1.3f\n" %(keht))
		
	fp.write("\ncell1\n")
	for x in range(natom_cell1):
			fp.write("%s %s %s %s\n"  %(ljust(str(cell1.atomType[x]),3), rjust(str(cell1.xyz[x][0]),20), rjust(str(cell1.xyz[x][1]),20), rjust(str(cell1.xyz[x][2]),20)) )
		

	fp.write("\ncell2\n")
	for x in range(natom_cell2):
			fp.write("%s %s %s %s\n"  %(ljust(str(cell2.atomType[x]),3), rjust(str(cell2.xyz[x][0]),20), rjust(str(cell2.xyz[x][1]),20), rjust(str(cell2.xyz[x][2]),20)) )
	
	fp.write("\n//////////////////////////////\n")
	fp.close

################################################
#
#	Export the husky file for the junction
#
################################################

def export_huckel(huckelFileName,cell1,cell2,hkl_param,keht):
	

	natom_cell1 = len(cell1.xyz)
	natom_cell2 = len(cell2.xyz)

	
	
	fp = open(huckelFileName,'w')
	fp.write(" //////////////////////////////\n //   pos for Huckel \n //////////////////////////////\n")
	fp.write("nb_atom %d\n" %(natom_cell1+natom_cell2) )
	fp.write("parameters %s\n" %(hkl_param))
	fp.write("Keht %1.3f\n" %(keht))
		
	
	for x in range(natom_cell1):
			fp.write("%s %s %s %s\n"  %(ljust(str(cell1.atomType[x]),3), rjust(str(cell1.xyz[x][0]),20), rjust(str(cell1.xyz[x][1]),20), rjust(str(cell1.xyz[x][2]),20)) )
		
	for x in range(natom_cell2):
			fp.write("%s %s %s %s\n"  %(ljust(str(cell2.atomType[x]),3), rjust(str(cell2.xyz[x][0]),20), rjust(str(cell2.xyz[x][1]),20), rjust(str(cell2.xyz[x][2]),20)) )
	
	fp.write("\n//////////////////////////////\n")
	fp.close


################################################
#
#	Export the eigenvector files
#
################################################

def write_vect(data,iorb,Ef,K,fname):

	# number of K points
	nKpoint = len(data)
	
	#  total number of orbitals
	norb = len(data[0])-1
	
	# number of band to be pritned
	nBand = len(iorb)
	
	# energy scale
	nE = 250
	E = np.linspace(-5,8,nE)
	
	# sigma for the Gaussian
	sigma=0.1
	
	
	f = open(fname,'w')
	gmax = 0
	
	# for all theKpoints
	for iK in range(nKpoint):
	
		# energy of the band
		e0 = data[iK][0] -Ef
		
		# for all the energy in the scale
		for iE in range(nE):
			c = 0
			
			# sum up all the coefficients +1 because [0] is the energy of the band
			for ip in range(len(iorb)):
				c += data[iK][iorb[ip]+1]
			
			# compute the Gaussian
			g = c * np.exp( -(E[iE]-e0)**2/(sigma*c) )
			
			if g>gmax:
				gmax = g
			
			# print in the file
			f.write('%f %f %f\n' %(K[iK],E[iE],g))
		f.write('\n')
	f.close()
	return gmax


################################################
#
#	Compute the charte dipole interaction
#
################################################
def charge_dipole(cell,vect,nb_orb):

	natom = len(cell.atomType)

	norb_tot = np.sum(nb_orb)
	Hint = np.zeros((norb_tot,norb_tot))

	kC = 1.44 # eV.Ang Coulomb constant
 
	Qpb = 2 	# charge on lead atoms
	Qi = -1 	# charge on iodide 
	D =  2 		# dipole moment
	eps =  2  	# dielectric constant

	# for all atoms in the cell
	for iA in range(natom):
		if cell.charge[iA] != 0:

			# position of the orbitals in the matrix
			if iA == 0:
				init_orb = 0
			else:
				init_orb = np.sum(nb_orb[:iA])+0

			norb_atom = nb_orb[iA]
			ind = range(init_orb,init_orb+4)

			for iD in range(cell.nb_dipoles):

				#print iD

				r = cell.xyz[iA,:]-cell.dipole_center[iD,:]
				dist = np.sqrt(np.sum((cell.xyz[iA,:]-cell.dipole_center[iD,:])**2))
				m = cell.dipole_orientation[iD,:]

				#r /= np.linalg.norm(r)
				#m /= np.linalg.norm(m)

				if cell.atomType[iA].lower() == 'pb':

					cd_int = -Qpb*D*np.dot(r,m)/(eps*dist**3)

					if 0:
						print " lead atom %d: dipole %d " %(iA, iD)
						print " \t center \t\t %f %f %f " %(cell.dipole_center[iD,0],cell.dipole_center[iD,1],cell.dipole_center[iD,2])
						print " \t orientation \t\t %f %f %f " %(cell.dipole_orientation[iD,0],cell.dipole_orientation[iD,1],cell.dipole_orientation[iD,2])
						print " \t interaction \t\t %f" %(cd_int)
						print ''
					

				if cell.atomType[iA].lower() == 'i':
					cd_int = -Qi*D*np.dot(r,m)/(eps*dist**3)

				if dist < 20:
					Hint[ind,ind] += cd_int

	return Hint

################################################
#
#	Compute the charte dipole interaction
#
################################################

def charge_dipole_v2(cell,vect,nb_orb):

	natom = len(cell.atomType)

	norb_tot = np.sum(nb_orb)
	Hint = np.zeros((norb_tot,norb_tot))

	kC = 1.44 # eV.Ang Coulomb constant
 
	Qpb = 2 	# charge on lead atoms
	Qi  = 2 	# charge on iodide 
	D =  2 		# dipole moment
	eps =  2  	# dielectric constant

	# for all atoms in the cell
	for iA in range(natom):
		if cell.charge[iA] != 0:

			# position of the orbitals in the matrix
			if iA == 0:
				init_orb = 0
			else:
				init_orb = np.sum(nb_orb[:iA])+0

			norb_atom = nb_orb[iA]
			ind = range(init_orb,init_orb+4)

			for iD in range(cell.nb_dipoles):

				#print iD

				r = cell.xyz[iA,:]-cell.dipole_center[iD,:]
				dist = np.sqrt(np.sum((cell.xyz[iA,:]-cell.dipole_center[iD,:])**2))
				m = cell.dipole_orientation[iD,:]

				#r /= np.linalg.norm(r)
				#m /= np.linalg.norm(m)

				if cell.atomType[iA].lower() == 'pb':

					cd_int = -Qpb*D*np.dot(r,m)/(eps*dist**3)

					if 0:
						print " lead atom %d: dipole %d " %(iA, iD)
						print " \t center \t\t %f %f %f " %(cell.dipole_center[iD,0],cell.dipole_center[iD,1],cell.dipole_center[iD,2])
						print " \t orientation \t\t %f %f %f " %(cell.dipole_orientation[iD,0],cell.dipole_orientation[iD,1],cell.dipole_orientation[iD,2])
						print " \t interaction \t\t %f" %(cd_int)
						print ''
					

				if cell.atomType[iA].lower() == 'i':
					cd_int = -Qi*D*np.dot(r,m)/(eps*dist**3)

				if dist < 20:
					Hint[ind,ind] += cd_int

	return Hint

################################################
#
#	Compute the effective mass at a given Kpt
#
################################################
def compute_effMass(kpt,eps_effMass,vect,hcell,scell,H_INT,S_INT,nelec):

	# determine the number of filled bands
	nband_filled = int(1.0*nelec/2.)
	E = []

	# size of the system
	norb = len(hcell)
	nCell = len(H_INT)

	# for all the kpoint we need
	# i.e. +/- eps and 2eps in each 3 direction
	# plus the kpt
	# compute the energies

	count_point = 0
	print_perc = 10

	nKpoints = 124

	for iX in range(-2,3):
		for iY in range(-2,3):
			for iZ in range(-2,3):



				if (1.0*count_point/(1.0*nKpoints) > 0.1):
					sys.stdout.write('\r\t\t\t\t\t'+str((print_perc)) + ' % done')
					sys.stdout.flush()
					count_point = 0
					print_perc += 10
					
				count_point += 1

				# form the current kpoint
				kcurr = kpt + eps_effMass*np.array([iX,iY,iZ])
				

				# value of the k-point in the brillouin zone
				kk = kcurr[0] * vect.rec1 + kcurr[1] * vect.rec2  + kcurr[2] * vect.rec3

				# Declare the Hamiltonian and include the reference cell
				H2D = np.zeros((norb,norb),dtype=complex)
				H2D += hcell
				
				# Declare the Overlap and include the reference cell
				S2D = np.zeros((norb,norb),dtype=complex)
				S2D += scell
				
				# include all the interaction we've found earlier
				for icell in range(nCell):
				
					# extract the data from the classes
					vcell = H_INT[icell].kpoint[0]*vect.vect1 + H_INT[icell].kpoint[1]*vect.vect2 + H_INT[icell].kpoint[2]*vect.vect3
					hint = H_INT[icell].matrix
					sint = S_INT[icell].matrix
					
					# add them to the Hamiltonian/overlap with the Bloch phase
					H2D += np.exp(1j*np.dot(kk,vcell))*hint
					S2D += np.exp(1j*np.dot(kk,vcell))*sint


				# compute the eigenvalues of the system
				alpha, beta, vl, vr, work, info = scila.lapack.zggev(H2D,S2D)
				w = alpha/beta

				# sort the eigenvalues
				ind_sort = np.argsort(w)
				w = np.sort(w)

				# store the valence and conduction energies
				E.append([[iX,iY,iZ],[w[nband_filled-1].real,w[nband_filled].real]])

	
	# form the  tensors 
	# of the second derivatives
	dE_val = np.zeros((3,3))
	dE_cond = np.zeros((3,3))

	

	# diagonal element of the valence tensor
	kdir = ['x','y','z']
	for idir in range(3):
		dE_val[idir,idir] = dEdx2(kdir[idir],E,eps_effMass,'val')
	

	# offdiagonal element of the valence tensor
	dE_val[0,1] = dEdxy('xy',E,eps_effMass,'val')
	dE_val[0,2] = dEdxy('xz',E,eps_effMass,'val')
	dE_val[1,2] = dEdxy('yz',E,eps_effMass,'val')
	dE_val[1,0] = dEdxy('yx',E,eps_effMass,'val')
	dE_val[2,0] = dEdxy('zx',E,eps_effMass,'val')
	dE_val[2,1] = dEdxy('zy',E,eps_effMass,'val')

	#print ''
	#print dE_val

	# diagonal element of the conduction tensor
	kdir = ['x','y','z']
	for idir in range(3):
		dE_cond[idir,idir] = dEdx2(kdir[idir],E,eps_effMass,'cond')

	# offdiagonal element of the conduction tensor
	dE_cond[0,1] = dEdxy('xy',E,eps_effMass,'cond')
	dE_cond[0,2] = dEdxy('xz',E,eps_effMass,'cond')
	dE_cond[1,2] = dEdxy('yz',E,eps_effMass,'cond')
	dE_cond[1,0] = dEdxy('yx',E,eps_effMass,'cond')
	dE_cond[2,0] = dEdxy('zx',E,eps_effMass,'cond')
	dE_cond[2,1] = dEdxy('zy',E,eps_effMass,'cond')

	# diagonalize the tensors
	Eval, Uval = scila.eigh(dE_val)
	Econd, Ucond = scila.eigh(dE_cond)

	# inverse of the eigenvalues
	mhole = 1./Eval
	melec = 1./Econd

	# return the eff mass only no eigenvector
	return mhole, melec


################################################
# auxiliary function for EffMass calcultaions 
#	Get the energy for specific Kpoint
################################################
def get_E_for_effMass(i,j,k,E,band):
	
	for iE in range(len(E)):
		if E[iE][0] == [i,j,k]:
			e = E[iE][1]
			break;
	if band == 'val':
		return e[0]

	elif band == 'cond':
		return e[1]

################################################
# auxiliary function for EffMass calcultaions 
#	Get value of the diagonal tensor element
# x = 'x','y','z' is the direction in kspace
################################################
def dEdx2(x,E,eps,band):
	dE = 0

	# if we want the x direction
	if x == 'x':
		dE -= (get_E_for_effMass(-2,0,0,E,band)+get_E_for_effMass(2,0,0,E,band))
		dE += 16 * (get_E_for_effMass(-1,0,0,E,band)+get_E_for_effMass(1,0,0,E,band))
		dE -= 30 * get_E_for_effMass(0,0,0,E,band)

	# if we want the y direction
	elif x == 'y':
		dE -= (get_E_for_effMass(0,-2,0,E,band)+get_E_for_effMass(0,2,0,E,band))
		dE += 16 * (get_E_for_effMass(0,-1,0,E,band)+get_E_for_effMass(0,1,0,E,band))
		dE -= 30 * get_E_for_effMass(0,0,0,E,band)

	# if we want the z direction
	elif x == 'z':
		dE -= (get_E_for_effMass(0,0,-2,E,band)+get_E_for_effMass(0,0,2,E,band))
		dE += 16 * (get_E_for_effMass(0,0,1,E,band)+get_E_for_effMass(0,0,1,E,band))
		dE -= 30 * get_E_for_effMass(0,0,0,E,band)

	return dE/12./eps/eps


################################################
# auxiliary function for EffMass calcultaions 
#	Get value of the offdiagonal tensor element
# xy = 'xy','yz','xz' ... are the direction in kspace
################################################
def dEdxy_h2(xy,E,eps,band):
	
	# template for the element of the numerical derivative
	index_xy = []
	index_xy.append([ -1, [ [1,0],[-1,0],[0,1],[0,-1] ] ])
	index_xy.append([  2, [ [0,0] ]  ])
	index_xy.append([  1, [ [1,1],[-1,-1] ] ])

	dE = 0
	if xy == 'xy':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(d1,d2,0,E,band)
		
	if xy == 'xz':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(d1,0,d2,E,band)

	if xy == 'yz':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(0,d1,d2,E,band)

	if xy == 'yx':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(d2,d1,0,E,band)
		
	if xy == 'zx':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(d2,0,d1,E,band)

	if xy == 'zy':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(0,d2,d1,E,band)

	return dE/2./eps/eps


################################################
# auxiliary function for EffMass calcultaions 
#	Get value of the offdiagonal tensor element
# xy = 'xy','yz','xz' ... are the direction in kspace
################################################
def dEdxy(xy,E,eps,band):
	
	# template for the element of the numerical derivative
	index_xy = []
	index_xy.append([ -63, [ [1,-2],[2,-1],[-2,1],[-1,2] ] ])
	index_xy.append([  63, [ [-1,-2],[-2,-1],[1,2],[2,1] ] ])
	index_xy.append([  44, [ [2,-2],[-2,2] ] ])
	index_xy.append([ -44, [ [-2,-2],[2,2] ] ])
	index_xy.append([  74, [ [-1,-1],[1,1] ] ])
	index_xy.append([ -74, [ [1,-1],[-1,1] ] ])

	dE = 0
	if xy == 'xy':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(d1,d2,0,E,band)
		
	if xy == 'xz':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(d1,0,d2,E,band)

	if xy == 'yz':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(0,d1,d2,E,band)

	if xy == 'yx':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(d2,d1,0,E,band)
		
	if xy == 'zx':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(d2,0,d1,E,band)

	if xy == 'zy':
		for iCoeff in range(len(index_xy)):
			for iDir in range(len(index_xy[iCoeff][1])):
				d1 = index_xy[iCoeff][1][iDir][0]
				d2 = index_xy[iCoeff][1][iDir][1]
				dE += index_xy[iCoeff][0]*get_E_for_effMass(0,d2,d1,E,band)

	return dE/600./eps/eps

























		
