//////////////////////
// header general 
////////////////////////

#ifndef read_input_hkl_gen_h
#define read_input_hkl_gen_h

void read_molecule_input_file(char *file_name, atom **MOL, int *no_atoms,
			    int **ind_cell_1, int *no_atm_cell1, int **ind_cell_2, int *no_atm_cell2, double *KEHT, char *PARAM);

void extract_atom_mol(atom *mol_only,  int *ind_atom_mol, int no_atom_mol, atom *MOL);
#endif