import numpy as np
import scipy.linalg as scila
import os
import sys
import argparse
from myclass import *
from tools import *
from compute_orbitals import *
import hkl_module as hkl
import ctypes


__diag__ = 'scipy'		# method to diagonalize the system
__print_vect__ = 0 		# do we print the projected band structure
__deb__ = 0 			# flag for debugging


IND_CELL_MAX = 3		# Maximum index from the reference cell
DIST_CUTOFF = 30.0		# distance cutoff for the itneractions

def main(argv):

	print "\n\n"
	print "========================================================" 
	print "==         pyBand 2.0 :                               ==" 
	print "==         Python Extended Huckel                     ==" 
	print "==         Band Structure Calculation                 ==" 
	print "======================================================== \n\n" 

  		
	##########################
	#
	# 	Parse the arguments
	#
	##########################
	
	parser = argparse.ArgumentParser(description='Compute the band structure of periodic crystals')
	
	#required (positional) arguments
	parser.add_argument('input', help = 'input file for the calculation')
	parser.add_argument('hkl_param', help = 'EHMO parameter for the calculations')
		
	#optional arguments
	parser.add_argument('-odir', '--output_directory', default='./',
											help = 'output directory where the file will be stored')


	parser.add_argument('-keht', '--KEHT', default=1.75,
											help = 'Keht constant for the calculation of the couplings',type=float)
		
	# done
	args=parser.parse_args()
	out_dir = args.output_directory
	
	##########################
	#
	# 	Read the data
	#
	##########################
	
	print " === Read the input file" 
	
	## read the position of the unit cell
	cell = unit_cell(args.input)
	
	## read the repetition vectors
	vect = vect_rep(args.input)
 	
	## read the bz path
	bzpath,dpath,dedge,kmesh = define_bzpath_vasp(args.input)

	## read the info for eff mass calculation
	__compute_effMass__, kpt_effMass, eps_effMass = read_effMass(args.input)
	nK_effMass = len(kpt_effMass)
	

	## determine the number of electrons and the # of orbtials per atoms in the cell
	nelec, nb_orb = count_electron(cell,args.hkl_param)
	

	## index of homo/lumo
	index_homo = nelec/2-1
	index_lumo = index_homo+1
	
	##########################
	#
	# 	Determine the space to
	#	  sample
	#
	##########################
	
	
	print " === Determine the K-points to sample"
	
	# default value for 3D system
	ind_min_1 = -IND_CELL_MAX
	ind_max_1 = IND_CELL_MAX

	ind_min_2 = -IND_CELL_MAX
	ind_max_2 = IND_CELL_MAX
	
	ind_min_3 = -IND_CELL_MAX
	ind_max_3 = IND_CELL_MAX
	
	# first direction must be non zero
	if sum(vect.vect1)==0:
		print " -- Warning the first repetition vector is null"
		
	# if the direction 2 is null (1D system)
	if sum(vect.vect2)<1E-16:
		ind_min_2 = 1
		ind_max_2 = 1
		print "Second vector null. Assuming a 1D system"

	# if the direction 3 is null (2D system)
	if sum(vect.vect3)<1E-16:
		ind_min_3 = 1
		ind_max_3 = 1
		print "Third vector null. Assuming a 2D system"

	##########################
	#
	# 	Compute the electronic
	#   Structure of the cells
	#   and their interactions
	#
	##########################
	
	# define the array of class that hold the interaction matrices
	# see myclass.py
	H_INT = []
	S_INT = []

	# variables for timing printer
	nb_cell_tot = (ind_max_1+1-ind_min_1)*(ind_max_1+1-ind_min_1)*(ind_max_3+1-ind_min_3)
	count_cell = 0
	print_perc = 10

	# read the numb of orbitals at the first iteration
	_read_orb_two_cells_ = 1 
	nb_orb_two_cells = 0
	nb_orb_one_cell = 0

	print " === Compute the electronic structure   ",

	# explore the 3d space
	for in1 in range(ind_min_1,ind_max_1+1):
		for in2 in range(ind_min_2,ind_max_2+1):
			for in3 in range(ind_min_3,ind_max_3+1):
				
				if (1.0*count_cell/(1.0*nb_cell_tot) > 0.1):
					#print " \r %d %% done" %(print_perc)
					sys.stdout.write('\r\t\t\t\t\t'+str(print_perc) + ' % done')
					sys.stdout.flush()
					count_cell = 0
					print_perc += 10
							
				# translation vector
				tr = in1*vect.vect1 + in2*vect.vect2 + in3*vect.vect3
			
				# distance to the reference unit cell
				dist_cell = np.sqrt(sum(tr**2))

				# if we're below the distance cutoff
				# but not on the reference cell
				if dist_cell < DIST_CUTOFF and dist_cell > 0.0:

					
					# copyt the cell and translate it
					cell2=copy.deepcopy(cell)
					cell2.translate(tr)

					# export the input file for huckel
					hkl_fn = "two_cell.in"
					export_huckel(hkl_fn,cell,cell2,args.hkl_param,args.KEHT)

					# compute the hamiltonian and overlap of the two cells
					if _read_orb_two_cells_:
						nb_orb_two_cells = hkl.nbOrb(hkl_fn)
						nb_orb_one_cell = nb_orb_two_cells/2
						_read_orb_two_cells_ = 0

					# initialize and compute the Hamiltonian and overlap
					h_twocells = np.zeros(nb_orb_two_cells*nb_orb_two_cells)
					s_twocells = np.zeros(nb_orb_two_cells*nb_orb_two_cells)
					hkl.hklHam(h_twocells,s_twocells,nb_orb_two_cells,hkl_fn)
				
					# reshape the matrices
					h_twocells = h_twocells.reshape(nb_orb_two_cells,nb_orb_two_cells)
					s_twocells = s_twocells.reshape(nb_orb_two_cells,nb_orb_two_cells)
					

					# extract the submatrices
					hcell = h_twocells[0:nb_orb_one_cell,0:nb_orb_one_cell]
					scell = s_twocells[0:nb_orb_one_cell,0:nb_orb_one_cell]
					

					hint = h_twocells[0:nb_orb_one_cell,nb_orb_one_cell:]
					sint = s_twocells[0:nb_orb_one_cell,nb_orb_one_cell:]

					# append it to the list of interaction matrices
					H_INT.append(int_mat(hint,[in1,in2,in3]))
					S_INT.append(int_mat(sint,[in1,in2,in3]))

					# increment the count
					count_cell += 1

	sys.stdout.write('\r\t\t\t\t\t'+str(100) + ' % done\n')
	os.system('rm two_cell.in')

	# number of cell in the system
	nCell = len(H_INT)

	# size pf the cell
	norb = len(hcell)

	# duplicate the dipoles to neighbouring cells
	#cell.duplicate_dipoles(vect)

	# include the charge-dipole interactions in hcell
	#print " === Compute the charge-dipole interactions   "
	#hchr_dip = charge_dipole_v2(cell,vect,nb_orb)
	#hcell += hchr_dip

	##############################################
	#
	#			Compute the Band 
	#			Structure
	#			Along the BZPATH
	#
	##############################################

	nKpoints = len(bzpath)
	W = zeros((nKpoints,norb))
	UH = zeros((nKpoints,norb+1))
	UL = zeros((nKpoints,norb+1))
	
	count_point = 0
	print_perc = 10
	print " === Compute the band structure",


	# walk along the path of the Brillouin zone
	for iKpoint in range(nKpoints):

	
		if (1.0*count_point/(1.0*nKpoints) > 0.1):
			sys.stdout.write('\r\t\t\t\t\t'+str(print_perc) + ' % done')
			sys.stdout.flush()
			count_point = 0
			print_perc += 10
			
		count_point += 1

		# value of the k-point in the brillouin zone
		kk = bzpath[iKpoint][0] * vect.rec1 + bzpath[iKpoint][1] * vect.rec2  + bzpath[iKpoint][2] * vect.rec3

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
			

		#################################################
		# diagonalize the generalized eigenvalue problem
		# with high level scipy routines 
		#################################################

		if __diag__ == 'scipy':
			# w are the eigenvalues
			# u are the left eigenvectors
			w,u = scila.eig(H2D,S2D)

		#################################################
		# diagonalize the generalized eigenvalue problem
		# with low level lapack routines 
		#################################################

		if __diag__ == 'lapack':
			# w are the eigenvalues
			# u are the left eigenvectors
			alpha, beta, vl, vr, work, info = scila.lapack.zggev(H2D,S2D)
			w = alpha/beta
			u = vl
		
		# sort the eigenvalues
		ind_sort = np.argsort(w)
		w = np.sort(w)

		# store the eigenvalues
		for iband in range(norb):
			W[iKpoint][iband] = real(w[iband])

		# sort the eigenvectors
		u = u[:,ind_sort]
		
		# Normalize the eigenvectors
		for i in range(len(u)):
			u[:,i] /= sqrt(dot(u[:, i], u[:, i]))
			#u[:,i] /= np.linalg.norm(u[:,i])

		# print the projection for debugging
		if __deb__ and iKpoint == 0:
			label = ['Pb:s','Pb:px','Pb:py','Pb:pz','I:s','I:px','I:py','I:pz','I:s','I:px','I:py','I:pz','I:s','I:px','I:py','I:pz']
			print "\n=== AOweight \t  HOMO \t\t  LUMO"
			for iX in range(len(u)):
				print '  %d %s  \t%1.2f %% \t\t%1.2f %%' %(iX,label[iX],abs(u[iX,index_homo])**2,abs(u[iX,index_lumo])**2)


		# form the matrix containing for each Kpoint
		# the energy and the AO coefficients of the HOMO and LUMO
		UH[iKpoint][0] = w[index_homo].real
		UL[iKpoint][0] = w[index_lumo].real
		for iorb in range(1,norb+1):
			UH[iKpoint][iorb] = abs(u[iorb-1,index_homo])**2
			UL[iKpoint][iorb] = abs(u[iorb-1,index_lumo])**2

		# compute and print the orbitals for 
		# each band at the edge of the bzpath
		if (iKpoint % kmesh) == 0 or iKpoint == nKpoints-1:
			filename = 'orbitals_k%d.dat' %(iKpoint)
			nb_print_orb = 20

			if nb_print_orb > norb:
				nb_print_orb = norb

			kpt_print = [bzpath[iKpoint][0],bzpath[iKpoint][1],bzpath[iKpoint][2]]
			print_mo_mopac(w.real,u.real,S2D.real,index_homo,cell,filename,args.hkl_param,nb_print_orb,kpt_print)
		
	sys.stdout.write('\r\t\t\t\t\t'+str(100) + ' % done\n')


	########################################
	##
	##		Compute the Effective
	##		Masses at the specified Kpoints
	##
	########################################
	print " === Compute the effective Mass   "

	if __compute_effMass__:
		EFFMASS_VAL = []
		EFFMASS_2eps_VAL = []
		EFFMASS_COND = []
		EFFMASS_2eps_COND = []
		for iKpoint in range(nK_effMass):

			print '   - %d/%d Kpoint (eps)' %(iKpoint+1,nK_effMass),
			effMass_val, effMass_cond = compute_effMass(kpt_effMass[iKpoint],eps_effMass,vect,hcell,scell,H_INT,S_INT,nelec)
			EFFMASS_VAL.append([kpt_effMass[iKpoint],effMass_val])
			EFFMASS_COND.append([kpt_effMass[iKpoint],effMass_cond])
			print ''

			print '   - %d/%d Kpoint (2 x eps)' %(iKpoint+1,nK_effMass),
			effMass_val_2eps, effMass_cond_2eps = compute_effMass(kpt_effMass[iKpoint],2*eps_effMass,vect,hcell,scell,H_INT,S_INT,nelec)
			EFFMASS_2eps_VAL.append([kpt_effMass[iKpoint],effMass_val_2eps])
			EFFMASS_2eps_COND.append([kpt_effMass[iKpoint],effMass_cond_2eps])
			print ''

			

	########################################
	##
	##		Print the results 
	##		in the output directory
	##
	########################################
	
	if (os.path.isdir(out_dir)):
		print " === Output directory %s already exist" %out_dir
	else:
		print " === Create output directory %s" %out_dir
		os.system('mkdir ' + out_dir)
	
	
	# determine the Fermi energy
	nband_filled = int(1.0*nelec/2.)
	if (nband_filled - 1.*nelec/2.) > 0.001:
		print " \t --- Warning non fully filled band. Number of electrons %d, number of bands %d" %(nelec,norb)
	Ef = max(W[:,nband_filled-1])
	
	
	# print the bands in a gnuplot file
	print " === Print the bands"
	f = open(out_dir + "/band_structure.dat",'w')
	for iband in range(norb):
		for iKpoint in range(nKpoints):
			f.write("%f %f\n" %(dpath[iKpoint],W[iKpoint][iband] - Ef))
		if iband>-1:
			f.write("\n")
	f.close()

	# print the band delimitation i.e. the symmetry points
	# in a gnuplot file
	f = open(out_dir + "/band_limit.dat",'w')
	for i in range(len(dedge)):
		f.write("%f -100\n" %(dedge[i]) )
		f.write("%f  100\n\n" %(dedge[i]) )
	f.close()

	# write the band edge etc
	f = open(out_dir + "/fermi_level.dat",'w')
	f.write("%f %f\n" %(min(dpath),Ef-Ef))
	f.write("%f %f\n\n" %(max(dpath),Ef-Ef))
	f.close()
	print " === Top valence band \t\t\t %1.3f eV (#%d)" %(Ef,nband_filled-1)
	print " === Bottom conduction band \t\t %1.3f eV (#%d)" %(min(W[:,nband_filled]),nband_filled)
	print " === Band gap \t\t\t\t % 1.3f eV" %(min(W[:,nband_filled])-Ef)


	# print the eigenvectors files
	if __print_vect__:
		ind_print = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
		ind_print_homo = [5,6,7,9,10,11,13,14,15]
		ind_print_homo = [0]
		ind_print_lumo = [5,6,7,9,10,11,13,14,15]
		hmax = write_vect(UH,ind_print_homo,Ef,dpath,out_dir +'/homo.dat')
		lmax = write_vect(UL,ind_print_lumo,Ef,dpath,out_dir +'/lumo.dat')
		print 'max values of HOMO %f \nmax value of LUMO %f\n' %(hmax,lmax)
		os.system("gnuplot plot_vect")

	# clean the files
	#os.system("rm H.dat S.dat hcell.dat hint.dat sint.dat scell.dat" );


	if __compute_effMass__:
		f = open(out_dir +'/effMass.dat','w')
		f.write(' Effective Mass Calculation \n\n')
		for iEff in range(len(EFFMASS_VAL)):

			f.write('   - k-Point coordinates (fractional) \t    %1.3f %1.3f %1.3f\n' %(EFFMASS_VAL[iEff][0][0],EFFMASS_VAL[iEff][0][1],EFFMASS_VAL[iEff][0][2]))
			f.write('   ------------------------------------------------------------------\n')
			f.write('     hole effective mass \n')
			f.write('	  eps = %f                             % 1.3f % 1.3f % 1.3f \n' %(eps_effMass, EFFMASS_COND[iEff][1][0],EFFMASS_COND[iEff][1][1],EFFMASS_COND[iEff][1][2]))
			f.write('	  eps = %f                             % 1.3f % 1.3f % 1.3f \n' %(2*eps_effMass, EFFMASS_2eps_COND[iEff][1][0],EFFMASS_2eps_COND[iEff][1][1],EFFMASS_2eps_COND[iEff][1][2]))
			f.write('   ------------------------------------------------------------------\n')
			f.write('     electron effective mass \n')
			f.write('	  eps = %f                             % 1.3f % 1.3f % 1.3f \n' %(eps_effMass, EFFMASS_VAL[iEff][1][0],EFFMASS_VAL[iEff][1][1],EFFMASS_VAL[iEff][1][2]))
			f.write('	  eps = %f                             % 1.3f % 1.3f % 1.3f \n' %(2*eps_effMass, EFFMASS_2eps_VAL[iEff][1][0],EFFMASS_2eps_VAL[iEff][1][1],EFFMASS_2eps_VAL[iEff][1][2]))
			f.write('\n\n')
		f.close()


	# plot the result
	os.chdir(out_dir)
	os.system("gnuplot plot")
	os.chdir('../')

	


	print "\n\n"
	print "========================================================" 
	print "==         Exit pyBand2.0_perovskite :                ==" 
	print "==         Calculation done with success              ==" 
	print "======================================================== \n\n" 

if __name__=='__main__':
	main(sys.argv[1:])




