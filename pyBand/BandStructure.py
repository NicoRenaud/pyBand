import numpy as np
import scipy.linalg as scila
import os
import sys
import argparse
import ctypes
from copy import deepcopy
import pickle

from .compute_orbitals import *
from . import hkl_module 
from .ElectronicSystem import *



class int_mat(object):

	'''
	interaction matrices for specific kpoints
	'''
	def __init__(self,m,k):
		self.matrix = m
		self.kpoint = k

class bzpath(object):

	'''
	hold the information about the bzpath
	'''
	def __init__(self,bzpath,dpath,dedge,kmesh):

		self.bzpath = bzpath
		self.dpath = dpath
		self.dedge = dedge
		self.kmesh = kmesh

class BandStructure(object):

	def __init__(self,input_file,hkl_param):


		print("\n\n")
		print("========================================================") 
		print("==         pyBand                                     ==") 
		print("==         Python Extended Huckel                     ==") 
		print("==         Band Structure Calculation                 ==") 
		print("======================================================== \n\n") 



		self.__diag__ = 'scipy'		# method to diagonalize the system
		self.__print_vect__ = 0 		# do we print the projected band structure
		self.__deb__ = 0 			# flag for debugging

		self.IND_CELL_MAX = 3		# Maximum index from the reference cell
		self.DIST_CUTOFF = 30.0		# distance cutoff for the itneractions

		self.keht = 1.875			# Huckel parameter

		self.input = input_file
		self.input_format = 'VASP'
		self.hkl_param = hkl_param
		
		## read the position of the unit cell
		self.cell = unit_cell(self.input)
		
		## read the repetition vectors
		self.vect = vect_rep(self.input)
	 	
		## read the bz path
		if self.input_format.upper() == 'VASP':
			self.bzpath = self._define_bzpath_vasp()
		elif self.input_format.upper() == 'ADF':
			self.bzpath = self._define_bzpath_adf()
		else:
			raise ValueError('Input format %s not supported' %self.input_format)

		## read the info for eff mass calculation
		self._read_effMass()
		self.nK_effMass = len(self.kpt_effMass)
		

		## determine the number of electrons and the # of orbtials per atoms in the cell
		self._count_electron()
		
		## index of homo/lumo
		self.index_homo = int(self.nelec/2-1)
		self.index_lumo = self.index_homo+1

		# determine the Kpoint boundary
		self._kpoints_boundary()

		self.out_dir = './'

	def _define_bzpath_adf(self):

		'''
		Read the bz path from the input file
		Weight the kpoint with distance (ADF)
		'''
		vect_raw = []
		
		# default value of kmesh
		kmesh = 3
		
		#open the input file, read the info, close the file
		inputFile1=open(self.input,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()

		#remove any newlines at the end of the input file
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1

		#make a matrix of the atom information 
		for x in range(lenInput):
			if inputFileLines[x].split("\n")[0]== "BZPATH":
				while inputFileLines[x+1].split("\n")[0] != "END":
					vect_raw.append(inputFileLines[x+1].split())
					x+=1
					
		for x in range(lenInput):
			if inputFileLines[x].split(" ")[0]== "KMESH":
				kmesh = int(inputFileLines[x].split(" ")[1])

		# declare the path
		nPoints = len(vect_raw)
		KPOINTS = zeros(( (nPoints-1)*kmesh,3))
		DIST = zeros(  (nPoints-1)*kmesh)
		DEDGE = zeros(nPoints)
		kk = 0
		d0 = 0
		nE = 1
		# define the path
		for iPoint in range(nPoints-1):
		
			x1 = linspace(float(vect_raw[iPoint][0]),float(vect_raw[iPoint+1][0]),kmesh)
			x2 = linspace(float(vect_raw[iPoint][1]),float(vect_raw[iPoint+1][1]),kmesh)
			x3 = linspace(float(vect_raw[iPoint][2]),float(vect_raw[iPoint+1][2]),kmesh)

			for iMesh in range(kmesh):

				KPOINTS[kk][0] = x1[iMesh]
				KPOINTS[kk][1] = x2[iMesh]
				KPOINTS[kk][2] = x3[iMesh]
				dist = sqrt(sum((KPOINTS[iPoint*kmesh]-KPOINTS[kk])**2))
				DIST[kk] = d0 + dist
				kk+=1
			d0 = DIST[kk-1]
			DEDGE[nE] = d0
			nE+=1

		return bzpath(KPOINTS,DIST,DEDGE,kmesh)


	def _define_bzpath_vasp(self):

		'''
		read the bz path frpm the input file
		do not weight the kpoints with distance (VASP)
		'''

		vect_raw = []
		
		# default value of kmesh
		kmesh = 3
		
		#open the input file, read the info, close the file
		inputFile1=open(self.input,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()

		#remove any newlines at the end of the input file
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1

		#make a matrix of the atom information 
		for x in range(lenInput):
			if inputFileLines[x].split("\n")[0]== "BZPATH":
				while inputFileLines[x+1].split("\n")[0] != "END":
					vect_raw.append(inputFileLines[x+1].split())
					x+=1
					
		for x in range(lenInput):
			if inputFileLines[x].split(" ")[0]== "KMESH":
				kmesh = int(inputFileLines[x].split(" ")[1])

		# declare the path
		nPoints = len(vect_raw)
		KPOINTS = zeros(( (nPoints-1)*kmesh,3))
		DIST = zeros(  (nPoints-1)*kmesh)
		DEDGE = zeros(nPoints)
		kk = 0
		d0 = 0
		nE = 1
		# define the path
		for iPoint in range(nPoints-1):
		
			x1 = linspace(float(vect_raw[iPoint][0]),float(vect_raw[iPoint+1][0]),kmesh)
			x2 = linspace(float(vect_raw[iPoint][1]),float(vect_raw[iPoint+1][1]),kmesh)
			x3 = linspace(float(vect_raw[iPoint][2]),float(vect_raw[iPoint+1][2]),kmesh)

			for iMesh in range(kmesh):

				KPOINTS[kk][0] = x1[iMesh]
				KPOINTS[kk][1] = x2[iMesh]
				KPOINTS[kk][2] = x3[iMesh]
				dist = sqrt(sum((KPOINTS[iPoint*kmesh]-KPOINTS[kk])**2))
				DIST[kk] = kk
				kk+=1
			d0 = DIST[kk-1]
			DEDGE[nE] = d0
			nE+=1

		return bzpath(KPOINTS,DIST,DEDGE,kmesh)

	def _read_effMass(self):

		'''
		read the information for eff mass calculation
		'''

		kpt_effMass = []
		__compute_effMass__ = 0
		eps_effMass = 0.0
		
		#open the input file, read the info, close the file
		inputFile1=open(self.input,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()

		#remove any newlines at the end of the input file
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1

		#look for the KPT key word and store the data
		for x in range(lenInput):
			if inputFileLines[x].split("\n")[0]== "KPT_EFFMASS":

				# change the flag and define the default value for eps
				__compute_effMass__ = 1
				eps_effMass = 0.01

				# read the kpt we want
				while inputFileLines[x+1].split("\n")[0] != "END":
					data = inputFileLines[x+1].split()
					kpt_effMass.append( [float(data[0]), float(data[1]), float(data[2]) ] )
					x+=1

		# find the eps value
		for x in range(lenInput):
			if inputFileLines[x].split(" ")[0]== "EPS_EFFMASS":
				eps_effMass = float(inputFileLines[x].split(" ")[1])


		self.kpt_effMass = kpt_effMass
		self.eps_effMass = eps_effMass


	def _count_electron(self):

		param_raw = []
		atom_type = []
		norb_type = []

		#open the input file, read the info, close the file
		inputFile1=open(self.hkl_param,'r')
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
			norb_type.append( [ int(param_raw[x][2]),  int(param_raw[x][3]),int(param_raw[x][4]), int(param_raw[x][5]) ] )

		nelec = 0
		nb_orb = []
		for i in range(self.cell.nAtom):
			for j in range(len(atom_type)):
				if self.cell.atomType[i].decode('utf-8') == atom_type[j][0]:

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
		self.nelec = nelec
		self.nb_orb = nb_orb

	def _kpoints_boundary(self):
	
		'''
		Determine the Kpoints to sample
		'''
		
		# default value for 3D system
		self.ind_min_1 = -self.IND_CELL_MAX
		self.ind_max_1 = self.IND_CELL_MAX

		self.ind_min_2 = -self.IND_CELL_MAX
		self.ind_max_2 = self.IND_CELL_MAX
		
		self.ind_min_3 = -self.IND_CELL_MAX
		self.ind_max_3 = self.IND_CELL_MAX
		
		# first direction must be non zero
		if sum(self.vect.vect1)==0:
			print(" -- Warning the first repetition vector is null")
			
		# if the direction 2 is null (1D system)
		if sum(self.vect.vect2)<1E-16:
			self.ind_min_2 = 1
			self.ind_max_2 = 1
			print("Second vector null. Assuming a 1D system")

		# if the direction 3 is null (2D system)
		if sum(self.vect.vect3)<1E-16:
			self.ind_min_3 = 1
			self.ind_max_3 = 1
			print("Third vector null. Assuming a 2D system")

	##################################################################################
	#
	#	Electronic Structure Calculation
	#
	##################################################################################

	def electronic_structure(self):

		'''
		Compute the electronic
		Structure of the cells
		and their interactions
		'''
		
		# define the array of class instance
		# for hamiltonian and Overlaperlap
		self.H_INT = []
		self.S_INT = []

		# variables for timing printer
		nb_cell_tot = (self.ind_max_1+1-self.ind_min_1)*(self.ind_max_1+1-self.ind_min_1)*(self.ind_max_3+1-self.ind_min_3)
		count_cell = 0
		print_perc = 10

		# read the numb of orbitals at the first iteration
		_read_orb_two_cells_ = 1 
		nb_orb_two_cells = 0
		nb_orb_one_cell = 0

		print(" === Compute the electronic structure   ", end=' ')

		# explore the 3d space
		for in1 in range(self.ind_min_1,self.ind_max_1+1):
			for in2 in range(self.ind_min_2,self.ind_max_2+1):
				for in3 in range(self.ind_min_3,self.ind_max_3+1):
					
					if (1.0*count_cell/(1.0*nb_cell_tot) > 0.1):
						#print " \r %d %% done" %(print_perc)
						sys.stdout.write('\r\t\t\t\t\t'+str(print_perc) + ' % done')
						sys.stdout.flush()
						count_cell = 0
						print_perc += 10
								
					# translation vector
					tr = in1*self.vect.vect1 + in2*self.vect.vect2 + in3*self.vect.vect3
				
					# distance to the reference unit cell
					dist_cell = np.sqrt(sum(tr**2))

					# if we're below the distance cutoff
					# but not on the reference cell
					if dist_cell < self.DIST_CUTOFF and dist_cell > 0.0:

						
						# copyt the cell and translate it
						cell2 = deepcopy(self.cell)
						cell2.translate(tr)

						# export the input file for huckel
						hkl_fn = b"two_cell.in"
						self._export_huckel(hkl_fn,self.cell,cell2)

						# compute the hamiltonian and overlap of the two cells
						if _read_orb_two_cells_:
							nb_orb_two_cells = hkl_module.nbOrb(hkl_fn)
							nb_orb_one_cell = int(nb_orb_two_cells/2)
							_read_orb_two_cells_ = 0

						# initialize and compute the Hamiltonian and overlap
						h_twocells = np.zeros(nb_orb_two_cells*nb_orb_two_cells)
						s_twocells = np.zeros(nb_orb_two_cells*nb_orb_two_cells)
						hkl_module.hklHam(h_twocells,s_twocells,nb_orb_two_cells,hkl_fn)
					
						# reshape the matrices
						h_twocells = h_twocells.reshape(nb_orb_two_cells,nb_orb_two_cells)
						s_twocells = s_twocells.reshape(nb_orb_two_cells,nb_orb_two_cells)
						

						# extract the submatrices
						self.hcell = h_twocells[0:nb_orb_one_cell,0:nb_orb_one_cell]
						self.scell = s_twocells[0:nb_orb_one_cell,0:nb_orb_one_cell]
						

						hint = h_twocells[0:nb_orb_one_cell,nb_orb_one_cell:]
						sint = s_twocells[0:nb_orb_one_cell,nb_orb_one_cell:]

						# append it to the list of interaction matrices
						self.H_INT.append(int_mat(hint,[in1,in2,in3]))
						self.S_INT.append(int_mat(sint,[in1,in2,in3]))

						# increment the count
						count_cell += 1

		sys.stdout.write('\r\t\t\t\t\t'+str(100) + ' % done\n')
		os.system('rm two_cell.in')

		# number of cell in the system
		self.nCell = len(self.H_INT)

		# size pf the cell
		self.norb = len(self.hcell)

		# duplicate the dipoles to neighbouring cells
		#cell.duplicate_dipoles(vect)

		# include the charge-dipole interactions in hcell
		#print " === Compute the charge-dipole interactions   "
		#hchr_dip = charge_dipole_v2(cell,vect,nb_orb)
		#hcell += hchr_dip


	def _export_huckel(self,huckelFileName,cell1,cell2):
		

		'''
		Export the husky file for the junction
		'''

		natom_cell1 = len(cell1.xyz)
		natom_cell2 = len(cell2.xyz)

		fp = open(huckelFileName,'w')
		fp.write(" //////////////////////////////\n //   pos for Huckel \n //////////////////////////////\n")
		fp.write("nb_atom %d\n" %(natom_cell1+natom_cell2) )
		fp.write("parameters %s\n" %(self.hkl_param))
		fp.write("Keht %1.3f\n" %(self.keht))
			
		
		for x in range(natom_cell1):
				fp.write("%s %s %s %s\n"  %( str(cell1.atomType[x].decode('utf-8')), str(cell1.xyz[x][0]), str(cell1.xyz[x][1]), str(cell1.xyz[x][2])))
			
		for x in range(natom_cell2):
				fp.write("%s %s %s %s\n"  %( str(cell2.atomType[x].decode('utf-8')), str(cell2.xyz[x][0]), str(cell2.xyz[x][1]), str(cell2.xyz[x][2])))
		
		fp.write("\n//////////////////////////////\n")
		fp.close

	##################################################################################
	#
	#	Compute the band structure
	#
	##################################################################################

	def compute_bands(self):

		'''
		Compute the band struture along the specified BZ path
		'''

		nKpoints = len(self.bzpath.bzpath)
		self.W = zeros((nKpoints,self.norb))
		self.UH = zeros((nKpoints,self.norb+1))
		self.UL = zeros((nKpoints,self.norb+1))
		
		count_point = 0
		print_perc = 10
		print(" === Compute the band structure", end=' ')


		# walk along the path of the Brillouin zone
		for iKpoint in range(nKpoints):

			if (1.0*count_point/(1.0*nKpoints) > 0.1):
				sys.stdout.write('\r\t\t\t\t\t'+str(print_perc) + ' % done')
				sys.stdout.flush()
				count_point = 0
				print_perc += 10
				
			count_point += 1

			# value of the k-point in the brillouin zone
			kk  = self.bzpath.bzpath[iKpoint][0] * self.vect.rec1 
			kk += self.bzpath.bzpath[iKpoint][1] * self.vect.rec2  
			kk += self.bzpath.bzpath[iKpoint][2] * self.vect.rec3

			# Declare the Hamiltonian and include the reference cell
			H2D = np.zeros((self.norb,self.norb),dtype=complex)
			H2D += self.hcell
			
			# Declare the Overlap and include the reference cell
			S2D = np.zeros((self.norb,self.norb),dtype=complex)
			S2D += self.scell
			
			# include all the interaction we've found earlier
			for icell in range(self.nCell):
			
				# extract the data from the classes
				vcell  = self.H_INT[icell].kpoint[0] * self.vect.vect1 
				vcell += self.H_INT[icell].kpoint[1] * self.vect.vect2 
				vcell += self.H_INT[icell].kpoint[2] * self.vect.vect3

				hint = self.H_INT[icell].matrix
				sint = self.S_INT[icell].matrix
				
				# add them to the Hamiltonian/overlap with the Bloch phase
				H2D += np.exp(1j*np.dot(kk,vcell))*hint
				S2D += np.exp(1j*np.dot(kk,vcell))*sint
				

			#################################################
			# diagonalize the generalized eigenvalue problem
			# with high level scipy routines 
			#################################################

			if self.__diag__ == 'scipy':
				# w are the eigenvalues
				# u are the left eigenvectors
				w,u = scila.eig(H2D,S2D)

			#################################################
			# diagonalize the generalized eigenvalue problem
			# with low level lapack routines 
			#################################################

			if self.__diag__ == 'lapack':
				# w are the eigenvalues
				# u are the left eigenvectors
				alpha, beta, vl, vr, work, info = scila.lapack.zggev(H2D,S2D)
				w = alpha/beta
				u = vl
			
			# sort the eigenvalues
			ind_sort = np.argsort(w)
			w = np.sort(w)

			# store the eigenvalues
			for iband in range(self.norb):
				self.W[iKpoint][iband] = real(w[iband])

			# sort the eigenvectors
			u = u[:,ind_sort]
			
			# Normalize the eigenvectors
			for i in range(len(u)):
				u[:,i] /= sqrt(dot(u[:, i], u[:, i]))
				#u[:,i] /= np.linalg.norm(u[:,i])

			# form the matrix containing for each Kpoint
			# the energy and the AO coefficients of the HOMO and LUMO
			self.UH[iKpoint][0] = w[self.index_homo].real
			self.UL[iKpoint][0] = w[self.index_lumo].real
			for iorb in range(1,self.norb+1):
				self.UH[iKpoint][iorb] = abs(u[iorb-1,self.index_homo])**2
				self.UL[iKpoint][iorb] = abs(u[iorb-1,self.index_lumo])**2

			'''
			# compute and print the orbitals for 
			# each band at the edge of the bzpath
			if (iKpoint % self.bzpath.kmesh) == 0 or iKpoint == nKpoints-1:
				filename = 'orbitals_k%d.dat' %(iKpoint)
				nb_print_orb = 20

				if nb_print_orb > self.norb:
					nb_print_orb = self.norb

				kpt_print = [self.bzpath.bzpath[iKpoint][0],self.bzpath.bzpath[iKpoint][1],self.bzpath.bzpath[iKpoint][2]]
				print_mo_mopac(w.real,u.real,S2D.real,self.index_homo,self.cell,filename,self.hkl_param,nb_print_orb,kpt_print)
			'''
			
		sys.stdout.write('\r\t\t\t\t\t'+str(100) + ' % done\n')

	##############################################################################
	#
	#	Compute the effecitve masses
	#
	###############################################################################
	def compute_effmass(self):

		'''
		Compute the effecitve masses at the desired Kpoints
		'''

		print(" === Compute the effective Mass   ")

		
		self.EFFMASS_VAL = []
		self.EFFMASS_2eps_VAL = []
		self.EFFMASS_COND = []
		self.EFFMASS_2eps_COND = []

		for iKpoint in range(self.nK_effMass):

			print('   - %d/%d Kpoint (eps)' %(iKpoint+1,self.nK_effMass), end=' ')
			effMass_val, effMass_cond = self._getmass(kpt_effMass[iKpoint],eps_effMass,vect,hcell,scell,H_INT,S_INT,nelec)
			self.EFFMASS_VAL.append([kpt_effMass[iKpoint],effMass_val])
			self.EFFMASS_COND.append([kpt_effMass[iKpoint],effMass_cond])
			print('')

			print('   - %d/%d Kpoint (2 x eps)' %(iKpoint+1,self.nK_effMass), end=' ')
			effMass_val_2eps, effMass_cond_2eps = self._getmass()
			self.EFFMASS_2eps_VAL.append([kpt_effMass[iKpoint],effMass_val_2eps])
			self.EFFMASS_2eps_COND.append([kpt_effMass[iKpoint],effMass_cond_2eps])
			print('')


	def export_effective_mass(self,fname):
		f = open(self.out_dir + fname,'w')
		f.write(' Effective Mass Calculation \n\n')
		for iEff in range(len(EFFMASS_VAL)):

			f.write('   - k-Point coordinates (fractional) \t    %1.3f %1.3f %1.3f\n' %(self.EFFMASS_VAL[iEff][0][0],self.EFFMASS_VAL[iEff][0][1],self.EFFMASS_VAL[iEff][0][2]))
			f.write('   ------------------------------------------------------------------\n')
			f.write('     hole effective mass \n')
			f.write('	  eps = %f                             % 1.3f % 1.3f % 1.3f \n' %(self.eps_effMass, self.EFFMASS_COND[iEff][1][0],self.EFFMASS_COND[iEff][1][1],self.EFFMASS_COND[iEff][1][2]))
			f.write('	  eps = %f                             % 1.3f % 1.3f % 1.3f \n' %(2*self.eps_effMass, self.EFFMASS_2eps_COND[iEff][1][0],self.EFFMASS_2eps_COND[iEff][1][1],self.EFFMASS_2eps_COND[iEff][1][2]))
			f.write('   ------------------------------------------------------------------\n')
			f.write('     electron effective mass \n')
			f.write('	  eps = %f                             % 1.3f % 1.3f % 1.3f \n' %(self.eps_effMass, self.EFFMASS_VAL[iEff][1][0],self.EFFMASS_VAL[iEff][1][1],self.EFFMASS_VAL[iEff][1][2]))
			f.write('	  eps = %f                             % 1.3f % 1.3f % 1.3f \n' %(2*self.eps_effMass, self.EFFMASS_2eps_VAL[iEff][1][0],self.EFFMASS_2eps_VAL[iEff][1][1],self.EFFMASS_2eps_VAL[iEff][1][2]))
			f.write('\n\n')
		f.close()


	def _getmass(self,iK):

		'''
		compute the effecitve mass for a Kpoint
		'''

		kpt = self.kpt_effMass[iK],

		# determine the number of filled bands
		nband_filled = int(1.0*self.nelec/2.)
		E = []

		# size of the system
		norb = len(self.hcell)
		nCell = len(self.H_INT)

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
					kcurr = kpt + self.eps_effMass*np.array([iX,iY,iZ])
					

					# value of the k-point in the brillouin zone
					kk = kcurr[0] * self.vect.rec1 + kcurr[1] * self.vect.rec2  + kcurr[2] * self.vect.rec3

					# Declare the Hamiltonian and include the reference cell
					H2D = np.zeros((self.norb,self.norb),dtype=complex)
					H2D += self.hcell
					
					# Declare the Overlap and include the reference cell
					S2D = np.zeros((self.norb,self.norb),dtype=complex)
					S2D += self.scell
					
					# include all the interaction we've found earlier
					for icell in range(nCell):
					
						# extract the data from the classes
						vcell = self.H_INT[icell].kpoint[0]*self.vect.vect1 + self.H_INT[icell].kpoint[1]*self.vect.vect2 + self.H_INT[icell].kpoint[2]*self.vect.vect3
						hint = self.H_INT[icell].matrix
						sint = self.S_INT[icell].matrix
						
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
			dE_val[idir,idir] = self._dEdx2(kdir[idir],E,eps_effMass,'val')
		

		# offdiagonal element of the valence tensor
		dE_val[0,1] = self._dEdxy('xy',E,eps_effMass,'val')
		dE_val[0,2] = self._dEdxy('xz',E,eps_effMass,'val')
		dE_val[1,2] = self._dEdxy('yz',E,eps_effMass,'val')
		dE_val[1,0] = self._dEdxy('yx',E,eps_effMass,'val')
		dE_val[2,0] = self._dEdxy('zx',E,eps_effMass,'val')
		dE_val[2,1] = self._dEdxy('zy',E,eps_effMass,'val')

		#print ''
		#print dE_val

		# diagonal element of the conduction tensor
		kdir = ['x','y','z']
		for idir in range(3):
			dE_cond[idir,idir] = self._dEdx2(kdir[idir],E,eps_effMass,'cond')

		# offdiagonal element of the conduction tensor
		dE_cond[0,1] = self._dEdxy('xy',E,eps_effMass,'cond')
		dE_cond[0,2] = self._dEdxy('xz',E,eps_effMass,'cond')
		dE_cond[1,2] = self._dEdxy('yz',E,eps_effMass,'cond')
		dE_cond[1,0] = self._dEdxy('yx',E,eps_effMass,'cond')
		dE_cond[2,0] = self._dEdxy('zx',E,eps_effMass,'cond')
		dE_cond[2,1] = self._dEdxy('zy',E,eps_effMass,'cond')

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
	#	Get value of the diagonal tensor element
	# x = 'x','y','z' is the direction in kspace
	################################################
	@staticmethod
	def _dEdx2(x,E,eps,band):
		dE = 0

		# if we want the x direction
		if x == 'x':
			dE -= (self._get_E_for_effMass(-2,0,0,E,band)+self._get_E_for_effMass(2,0,0,E,band))
			dE += 16 * (self._get_E_for_effMass(-1,0,0,E,band)+self._get_E_for_effMass(1,0,0,E,band))
			dE -= 30 * self._get_E_for_effMass(0,0,0,E,band)

		# if we want the y direction
		elif x == 'y':
			dE -= (self._get_E_for_effMass(0,-2,0,E,band)+self._get_E_for_effMass(0,2,0,E,band))
			dE += 16 * (self._get_E_for_effMass(0,-1,0,E,band)+self._get_E_for_effMass(0,1,0,E,band))
			dE -= 30 * self._get_E_for_effMass(0,0,0,E,band)

		# if we want the z direction
		elif x == 'z':
			dE -= (self._get_E_for_effMass(0,0,-2,E,band)+self._get_E_for_effMass(0,0,2,E,band))
			dE += 16 * (self._get_E_for_effMass(0,0,1,E,band)+self._get_E_for_effMass(0,0,1,E,band))
			dE -= 30 * self._get_E_for_effMass(0,0,0,E,band)

		return dE/12./eps/eps


	################################################
	# auxiliary function for EffMass calcultaions 
	#	Get value of the offdiagonal tensor element
	# xy = 'xy','yz','xz' ... are the direction in kspace
	################################################
	@staticmethod
	def _dEdxy_h2(xy,E,eps,band):
		
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
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(d1,d2,0,E,band)
			
		if xy == 'xz':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(d1,0,d2,E,band)

		if xy == 'yz':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(0,d1,d2,E,band)

		if xy == 'yx':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(d2,d1,0,E,band)
			
		if xy == 'zx':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(d2,0,d1,E,band)

		if xy == 'zy':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(0,d2,d1,E,band)

		return dE/2./eps/eps


	################################################
	# auxiliary function for EffMass calcultaions 
	#	Get value of the offdiagonal tensor element
	# xy = 'xy','yz','xz' ... are the direction in kspace
	################################################
	@staticmethod
	def _dEdxy(xy,E,eps,band):
		
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
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(d1,d2,0,E,band)
			
		if xy == 'xz':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(d1,0,d2,E,band)

		if xy == 'yz':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(0,d1,d2,E,band)

		if xy == 'yx':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(d2,d1,0,E,band)
			
		if xy == 'zx':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(d2,0,d1,E,band)

		if xy == 'zy':
			for iCoeff in range(len(index_xy)):
				for iDir in range(len(index_xy[iCoeff][1])):
					d1 = index_xy[iCoeff][1][iDir][0]
					d2 = index_xy[iCoeff][1][iDir][1]
					dE += index_xy[iCoeff][0]*self._get_E_for_effMass(0,d2,d1,E,band)

		return dE/600./eps/eps

	################################################
	# auxiliary function for EffMass calcultaions 
	#	Get the energy for specific Kpoint
	################################################
	@staticmethod
	def _get_E_for_effMass(i,j,k,E,band):
		
		for iE in range(len(E)):
			if E[iE][0] == [i,j,k]:
				e = E[iE][1]
				break;
		if band == 'val':
			return e[0]

		elif band == 'cond':
			return e[1]


	############################################################################
	def plot_bands(self,fermi=True,minE=None,maxE=None):

		import matplotlib.pyplot as plt 
		
		nband_filled = int(1.0*self.nelec/2.)
		Ef = max(self.W[:,nband_filled-1])

		if fermi:
			self.W-=Ef
			Ef = 0

		if minE is None:
			minE = np.min(self.W)-1.0

		if maxE is None:
			maxE = np.max(self.W)+1.0

		# plot occupied bands
		for iband in range(nband_filled):
			plt.plot(self.bzpath.dpath,self.W[:,iband],color='blue')

		# plot virtual bands
		for iband in range(nband_filled,self.norb):
			plt.plot(self.bzpath.dpath,self.W[:,iband],color='red')

		# plot edges
		for i in range(len(self.bzpath.dedge)):
			plt.plot([self.bzpath.dedge[i]]*2, [minE,maxE],color='black')


		# plot Ef
		plt.plot(self.bzpath.dpath[[0,-1]],[Ef,Ef],color='black')

		plt.ylim((minE,maxE))
		plt.show()

	def pickle(self,fname='./bands.pkl'):
		import pickle
		pickle.dump(self,open(fname,"wb"))






