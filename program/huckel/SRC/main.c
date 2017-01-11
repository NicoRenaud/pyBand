#include "./defMacro.h"
#include "./algebra_real/algebra.h"
#include "./algebra_real/extract_hamiltonian.h"
#include "./algebra_cplx/algebra_cplx.h"
#include "./huckel/header.h"
#include "./huckel/huckel.h"
#include "./huckel/read_input_hkl_general.h"


///////////////////////////////////////////////////////////
////////                                        ///////////
////////	MAIN FUNCTION                         ///////////
////////                                        ///////////
///////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  
  ////////////////////////
  //	DECLARATION	//
  ////////////////////////
  
  // data of the computation
  //SYS_XYZ  sys;
  
  
  // molecular Hamiltonian
  double *H,*S;

  double *H1,*S1;
  double *H12,*S12;
  
  // cluster/molecule/cluster information
  atom *SYSTEM;
  int nb_atom_tot;
  int nb_orb = 0;
  
  
  // molecule information
  atom *CELL1;

  
  // first cell
  int nb_atom_cell1 = 0;
  int nb_orb_cell_1 = 0;
  int *index_orb_cell_1;
  int *index_cell_1;
  

  // second cell
  int nb_atom_cell2;
  int nb_orb_cell_2 = 0;
  int *index_orb_cell_2;
  int *index_cell_2;

  // out files
  char  out_file_system[200];
  strcpy(out_file_system,argv[2]);
  strcat(out_file_system,"system.xyz");
	
	// out file for the matrix
	char out_hcell[200];
	strcpy(out_hcell,argv[2]);
  strcat(out_hcell,"hcell.dat");
	
	// out file for the matrix
	char out_scell[200];
	strcpy(out_scell,argv[2]);
  strcat(out_scell,"scell.dat");
	
	// out file for the matrix
	char out_hint[200];
	strcpy(out_hint,argv[2]);
  strcat(out_hint,"hint.dat");
	
	// out file for the matrix
	char out_sint[200];
	strcpy(out_sint,argv[2]);
  strcat(out_sint,"sint.dat");
  
  
  ////////////////////////////////////////////////////////////////////
  // Compute the total Hamiltonian in extended Huckel
  ////////////////////////////////////////////////////////////////////
     
  
  // compute the total Hamiltonian
  compute_huckel_hamiltonian_general(&H, &S, &nb_orb, &nb_atom_tot, &SYSTEM,
				&index_cell_1, &nb_atom_cell1,  &index_cell_2, &nb_atom_cell2,
				&index_orb_cell_1,&nb_orb_cell_1,&index_orb_cell_2, &nb_orb_cell_2,
				argv[1], argv[2]);

	
  ////////////////////////////////////////////////////////////////////
  // Extract the Hamiltonian of each fragments
  ////////////////////////////////////////////////////////////////////
  
  // extract the Hamiltonian/Overlap of the first cell
  H1 = calloc(nb_orb_cell_1*nb_orb_cell_1,sizeof(double));
  extract(H1, H, index_orb_cell_1, nb_orb_cell_1, index_orb_cell_1,nb_orb_cell_1, nb_orb, out_hcell);
  
  S1 = calloc(nb_orb_cell_1*nb_orb_cell_1,sizeof(double));
  extract(S1, S, index_orb_cell_1, nb_orb_cell_1, index_orb_cell_1,nb_orb_cell_1, nb_orb,out_scell);
	
	// extract the interaction Hamiltonian/Overlap
  H12 = calloc(nb_orb_cell_1*nb_orb_cell_2,sizeof(double));
  extract(H12, H, index_orb_cell_1, nb_orb_cell_1, index_orb_cell_2,nb_orb_cell_2, nb_orb,out_hint);

  S12 = calloc(nb_orb_cell_1*nb_orb_cell_2,sizeof(double));
  extract(S12, S, index_orb_cell_1, nb_orb_cell_2, index_orb_cell_2,nb_orb_cell_2, nb_orb,out_sint);


  ///////////////////////////////////
  //	free memory
  ///////////////////////////////////
  free(H1);
  free(S1);
  free(H12);
  free(S12);
  free(CELL1);
  
  
   ////////////////////////////////////////////////
  // 	   	Say Good bye			//
  /////////////////////////////////////////////// 
  //printf("\n====  ===============================  ====\n");
  //printf("====  Thanks for using Huckel 	       ====\n");
  //printf("====                    see you later  ====\n");
  //printf("====  ===============================  ====\n");

  return(0);
}
 
 
