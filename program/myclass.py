from numpy import *
import sys

##############################################
# class of the repetition vector
# and dimensionality of the system
##############################################

class vect_rep:

	# define the vector from the inut file
	def __init__(self,file_name):
	
		# define the raw vectors
		vect_raw = []
	
		#open the input file, read the info, close the file
		inputFile1=open(file_name,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()


		#remove any newlines at the end of the input file
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1

		#make a matrix of the atom information (this happens quite a bit in this code - get used to it)
		for x in range(lenInput):
			if inputFileLines[x].split("\n")[0]== "LATTICE":
				while inputFileLines[x+1].split("\n")[0] != "END":
					vect_raw.append(inputFileLines[x+1].split())
					x+=1
					
		# determine the size of the system
		nDim = len(vect_raw)

		# store the dimensionality
		vect_rep.dim = nDim
		
		# declare the vectors
		vect_rep.vect1 = zeros(3)
		vect_rep.vect2 = zeros(3)
		vect_rep.vect3 = zeros(3)
		
		# fill them in
		if nDim>=1:
			for iSub in range(nDim):
				vect_rep.vect1[iSub] = float(vect_raw[0][iSub])
		if nDim>=2:
			for iSub in range(nDim):
				vect_rep.vect2[iSub] = float(vect_raw[1][iSub])
		if nDim>=3:
			for iSub in range(nDim):
				vect_rep.vect3[iSub] = float(vect_raw[2][iSub])

		# declare the reciprocal lattice
		vect_rep.rec1 = zeros(3)
		vect_rep.rec2 = zeros(3)
		vect_rep.rec3 = zeros(3)

		if nDim==3:
			vect_rep.rec1 = 2.*pi*cross(vect_rep.vect2,vect_rep.vect3)/(dot(vect_rep.vect1,cross(vect_rep.vect2,vect_rep.vect3)))
			vect_rep.rec2 = 2.*pi*cross(vect_rep.vect3,vect_rep.vect1)/(dot(vect_rep.vect1,cross(vect_rep.vect2,vect_rep.vect3)))
			vect_rep.rec3 = 2.*pi*cross(vect_rep.vect1,vect_rep.vect2)/(dot(vect_rep.vect1,cross(vect_rep.vect2,vect_rep.vect3)))

		if nDim==2:
			vect_rep.vect3[2] = -1
			vect_rep.rec1 = 2.*pi*cross(vect_rep.vect2,vect_rep.vect3)/(dot(vect_rep.vect1,cross(vect_rep.vect2,vect_rep.vect3)))
			vect_rep.rec2 = 2.*pi*cross(vect_rep.vect3,vect_rep.vect1)/(dot(vect_rep.vect1,cross(vect_rep.vect2,vect_rep.vect3)))
			vect_rep.vect3[2] = 0

		if nDim==1:
			vect_rep.rec1 = 2.*pi*vect_rep.vect1
			

##############################################
# class of the repetition vector
# and dimensionality of the system
##############################################
class unit_cell:

	
	# define the molecule from the input files
	# file_name contaons the xyz

	def __init__(self,file_name):

		####################################################
		##
		## 	Read the XYZ Information
		##
		####################################################

		vect_raw = []

		#open the input file, read the info, close the file
		inputFile1=open(file_name,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()


		#remove any newlines at the end of the input file
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1

		#make a matrix of the atom information (this happens quite a bit in this code - get used to it)
		for x in range(lenInput):
			if inputFileLines[x].split("\n")[0]== "ATOMS":
				while inputFileLines[x+1].split("\n")[0] != "END":
					vect_raw.append(inputFileLines[x+1].split())
					x+=1

		# number of atoms
		self.nAtom = int(len(vect_raw))

		
		# type of atoms
		self.atomType = chararray(self.nAtom,itemsize=2)
		self.nb_dipoles = 0
		self.charge = zeros(self.nAtom)
		indexC = []
		indexN = []
		for i in range(self.nAtom):
			self.atomType[i] = vect_raw[i][0]
			if self.atomType[i].lower() == 'c':
				self.nb_dipoles += 1
				indexC.append(i)
			if self.atomType[i].lower() == 'n':
				indexN.append(i)
			if self.atomType[i].lower() == 'pb':
				self.charge[i] = +2
			if self.atomType[i].lower() == 'i':
				self.charge[i] = -1

		# atom positions	
		self.xyz = zeros((self.nAtom,3))	
		for i in range(self.nAtom):
			pos = [ float(vect_raw[i][1]) , float(vect_raw[i][2]),float(vect_raw[i][3])  ]
			self.xyz[i] = pos

		# dipoles
		'''
		self.dipole_center = zeros((self.nb_dipoles,3))
		self.dipole_orientation = zeros((self.nb_dipoles,3))
		for i in range(len(indexC)):
			self.dipole_center[i,:] = 0.5*(self.xyz[indexC[i]]+self.xyz[indexN[i]])
			self.dipole_orientation[i,:] = self.xyz[indexC[i]]-self.xyz[indexN[i]]
		'''
	##############################
	#
	# translate the molecule
	#
	##############################
	def translate(self,tr):

		for i in range(self.nAtom):
			self.xyz[i] += tr

	##############################
	#
	# translate the molecule
	#
	##############################
	def duplicate_dipoles(self,vect):

		#nb_supp = 26*self.nb_dipoles
		new_dipoles_center = [] #zeros((nb_supp,3))
		new_dipoles_orientation = [] #zeros((nb_supp,3))
		N = 5
		RANGE = range(-N,N+1)
		for ix in RANGE:
			for iy in RANGE:
				for iz in RANGE:
					if ix != 0 or iy != 0 or iz !=0:
						tr = ix * vect.vect1 + iy*vect.vect2 + iz*vect.vect3
						for iD in range(self.nb_dipoles):
							new_dipoles_center.append(self.dipole_center[iD,:] + tr)
							new_dipoles_orientation.append(self.dipole_orientation[iD,:])
							
		self.dipole_center = concatenate((self.dipole_center,new_dipoles_center))
		self.dipole_orientation = concatenate((self.dipole_orientation,new_dipoles_orientation))
		self.nb_dipoles = len(self.dipole_center)
		

##############################################
#  INTERACTIONS MATRICES
##############################################
class int_mat(object):

	def __init__(self,m,k):
		self.matrix = m
		self.kpoint = k

##############################################
# read the bz path frpm the input file
# weight the kpoint with distance (ADF)
##############################################
def define_bzpath_adf(file_name):

		vect_raw = []
		
		# default value of kmesh
		kmesh = 3
		
		#open the input file, read the info, close the file
		inputFile1=open(file_name,'r')
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

		return KPOINTS,DIST,DEDGE,kmesh

##############################################
# read the bz path frpm the input file
# do not weight the kpoints with distance (VASP)
##############################################
def define_bzpath_vasp(file_name):

		vect_raw = []
		
		# default value of kmesh
		kmesh = 3
		
		#open the input file, read the info, close the file
		inputFile1=open(file_name,'r')
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

		return KPOINTS,DIST,DEDGE,kmesh

##############################################
# read the information for eff mass calculation
##############################################
def read_effMass(file_name):

	kpt_effMass = []
	__compute_effMass__ = 0
	eps_effMass = 0.0
	
	#open the input file, read the info, close the file
	inputFile1=open(file_name,'r')
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


	return __compute_effMass__, kpt_effMass, eps_effMass
