#include "../defMacro.h"
#include "./header.h"
#include "./read_input_hkl.h"


void read_molecule_input_file(char *file_name, atom **MOL, int *no_atoms,
			    int **ind_cell_1, int *no_atm_cell1, int **ind_cell_2, int *no_atm_cell2, double *KEHT, char *PARAM )
{
  
  
  int i;
  int INT;
  int INT2;
  float X,Y,Z;
  
  int current_cell = 0;
  int read_cell[2];
  
  int no_atom_tot=-1;
  int nb_atom_tot_read=0;
  
		
  int *no_atom_cell;
  no_atom_cell = calloc(2,sizeof(int));
  no_atom_cell[0]  = -1;
  no_atom_cell[1]  = -1;
  
  int *no_atom_cell_read;
  no_atom_cell_read = calloc(2,sizeof(int));
  no_atom_cell_read[0]  = 0;
  no_atom_cell_read[1]  = 0;
  

  
  
  atom *MOL_TEMP=NULL;
  atom *MORE_MOL;
  
  char str_cmt[5];
  char str[100];
	

  int *index_cell_1=NULL;
  int *index_cell_2=NULL;
  
  int *more_index_cell_1;
  int *more_index_cell_2;
  
  int index_check;
  
  FILE *f;
  
  
  // check if the file exists
   f = fopen(file_name,"r");
   if (!f) {
      printf("couldn't read input file name %s\n", file_name);
      exit(1);
   }
   
   //////////////////////////////////////////
   // for all lines of the input file
   //////////////////////////////////////////
   while (!feof(f))
   //while(fgets(str,20,f) != NULL)
   {
     
        ///////////////////////
	// read the string
	///////////////////////
	fscanf(f,"%s",str);
	

	
	////////////////////////
	// if it's a comment
	////////////////////////
	strncpy(str_cmt,str,2);
	str_cmt[2] = '\0';
	if (!strcmp(str_cmt,"//"))
	  fgets(str,100,f);
	
	///////////////////////////
	// if it's not a comment
	///////////////////////////
	else
	{
	  
 ////////////////////////////////////////////////////////////////////////////////
 //				GENERAL INFORMATIONS
 ////////////////////////////////////////////////////////////////////////////////
 
	    ////////////////////////
	    // total number of atoms
	    ///////////////////////    
	    if( !strcmp(str,"nb_atom_tot") || !strcmp(str,"NB_ATOM_TOT")   || 
		  !strcmp(str,"nb_atom")  || !strcmp(str,"NB_ATOM")  ) 
	    {
        fscanf(f,"%d",&no_atom_tot);		  
        MOL_TEMP = calloc(no_atom_tot,sizeof(atom));
	    }  
		  
	    /////////////////////////////////
	    // number of atoms in the cell
      /////////////////////////////////
	    else  if( !strcmp(str,"nb_atom_cell") )
      {
        fscanf(f,"%d %d",&INT,&INT2);

          //no_atom_elec[INT-1]=INT2;
        
          if(INT == 1)
          {
            no_atom_cell[0]=INT2;
            index_cell_1 = calloc(no_atom_cell[0],sizeof(int));

          }
          else if(INT == 2)
          {
            no_atom_cell[1]=INT2;
            index_cell_2 = calloc(no_atom_cell[1],sizeof(int));
          }

	     }
			 
	    ////////////////////////
	    // mixing term for ehmo
	    ///////////////////////
	    else if( !strcmp(str,"Keht") || !strcmp(str,"KEHT")   ||
		  !strcmp(str,"keht")  || !strcmp(str,"kEHT")  ) 
	    {
        fscanf(f,"%lf",KEHT);
	    }
			
	    ////////////////////////
	    // EHMO parameters
			///////////////////////    
	    else if( !strcmp(str,"parameters") || !strcmp(str,"PARAMETERS")   ||
		  !strcmp(str,"param")  || !strcmp(str,"PARAM")  ) 
	    {
        fscanf(f,"%s",PARAM);
	    }
	     
      
            
	     
 ////////////////////////////////////////////////////////////////////////////////
 //				WHICH POSITIONS ARE WE READING
 ///////////////////////////////////////////////////////////////////////////////	
		  
	    ////////////////////////
	    // electrode
	    ///////////////////////    
	     else if( !strcmp(str,"cell1") || !strcmp(str,"cell2")  )
	     {
         // read the current elec
         //fscanf(f,"%d",&current_elec);
         current_cell++;
         
         // reinit the value of read_X
         for(i=0;i<2;i++)
           read_cell[i] = 0;
		  
         // init value of read_elec(current)
         read_cell[current_cell-1]  = 1;
         
	     }
	     	     
 ////////////////////////////////////////////////////////////////////////////////
 //			READING POSITIONS
 ///////////////////////////////////////////////////////////////////////////////
            
	    else
	    {
				 
	       index_check = findIndex_reduced(str,PARAM);
				 //printf("str : %s index : %d\n",str,index_check);
	       if (index_check > 0)
	       {
		
					// read the data
					fscanf(f,"%f%f%f",&X,&Y,&Z);
					
					// realloc the memory if not dine previously
					if(no_atom_tot == -1 || nb_atom_tot_read >= no_atom_tot)
					{
						MORE_MOL = realloc(MOL_TEMP,(nb_atom_tot_read+1)*sizeof(atom));
						MOL_TEMP = MORE_MOL;
					}
					
					// store the data
					strcpy(MOL_TEMP[nb_atom_tot_read].atomTypeChar,str);
					MOL_TEMP[nb_atom_tot_read].atomtype = index_check;
					MOL_TEMP[nb_atom_tot_read].x=X;
					MOL_TEMP[nb_atom_tot_read].y=Y;
					MOL_TEMP[nb_atom_tot_read].z=Z;
					
					
					// update the  atoms and their indexes of the first cell
					if(read_cell[0])
					{
						
						// realloc the indexes if not done previously
						if(no_atom_cell[0]==-1 || no_atom_cell_read[0]>= no_atom_cell[0])
						{
							more_index_cell_1 = realloc( index_cell_1,(no_atom_cell_read[0]+1)*sizeof(int));
							index_cell_1 = more_index_cell_1;
						}
						
						// store the index
							index_cell_1[ no_atom_cell_read[0] ] = nb_atom_tot_read;
							
							
						// one more atom in the first elec
						no_atom_cell_read[0] ++;
					}
					
					// update the  atoms and their indexes of the second cell
					if(read_cell[1])
					{
						
						// realloc the indexes if not done previously
						if(no_atom_cell[1]==-1 || no_atom_cell_read[1]>= no_atom_cell[1])
						{
							more_index_cell_2 = realloc(index_cell_2,(no_atom_cell_read[1]+1)*sizeof(int));
							index_cell_2 = more_index_cell_2;
						}
						
						// store the index
							index_cell_2[ no_atom_cell_read[1] ] = nb_atom_tot_read;
							
							
						// one more atom in the secon elec
						no_atom_cell_read[1] ++;
					}
						
					// one more atom read!
					nb_atom_tot_read ++;
		
		 }
		}
	 }
	}

      // export the data
      *MOL = MOL_TEMP;
      *ind_cell_1 = index_cell_1;
      *ind_cell_2 = index_cell_2;
      *no_atoms = nb_atom_tot_read;
		  //*KEHT = keht_temp;
      *no_atm_cell1 = no_atom_cell_read[0];
      *no_atm_cell2 = no_atom_cell_read[1];
      
      if(nb_atom_tot_read != no_atom_tot)
	printf("\nWarning: input claims there is %d atom in the system but %d have been found \n Please check %s\n",no_atom_tot,nb_atom_tot_read,file_name);
         

      if( no_atom_cell_read[0]!= no_atom_cell[0])
	printf("\nWarning: input claims there is %d atoms in the cell 1 but %d have been found \n Please check %s\n",no_atom_cell[0],no_atom_cell_read[0],file_name);
      
       if( no_atom_cell_read[1]!= no_atom_cell[1])
	printf("\nWarning: input claims there is %d atoms in cell 2 but %d have been found \n Please check %s\n",no_atom_cell[1],no_atom_cell_read[1],file_name);
       
      // close the file
      fclose(f);

}




void extract_atom_mol(atom *mol_only,  int *ind_atom_mol, int no_atom_mol, atom *MOL)
{
  
  int i;
  for(i=0;i<no_atom_mol;i++)
  {
   	strcpy(mol_only[i].atomTypeChar,MOL[ind_atom_mol[i]].atomTypeChar);
	mol_only[i].atomtype = MOL[ind_atom_mol[i]].atomtype;
	mol_only[i].x = MOL[ind_atom_mol[i]].x;
	mol_only[i].y = MOL[ind_atom_mol[i]].y;
	mol_only[i].z = MOL[ind_atom_mol[i]].z; 
  }
  
}
