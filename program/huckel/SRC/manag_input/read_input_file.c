#include "../defMacro.h"


/////////////////////////////////////////////////////////////
//	read the input file and calls
//	and generate the data for a junction
/////////////////////////////////////////////////////////////
void read_data(SYS_XYZ* sys_xyz, char *file_name)
{

   int i;
   int INT;
   char str [200];
   char str2[200];
   char str_out[200];
   char str_cmt[3];
   FILE *f;

   f = fopen(file_name,"r");
   if (!f) {
      printf("couldn't read %s\n", file_name);
      exit(1);
   }

   ///////////////////////////
   // default value
   ///////////////////////////
   
  
   // use which orbitals for TE all s_only, p_only px_only ...
   strcpy(sys_xyz->orb_contact,"all");
   
   // type of electronic structure
   strcpy(sys_xyz->elec_struct,"ehmo");
   
   // use all the orbital in the electrode all, s_only ....
   strcpy(sys_xyz->orb_elec,"all");
  
  // expor the MO
   strcpy(sys_xyz->export_MO,"yes");
   
   


   
   //////////////////////////////////////////
   // for all lines of the input file
   //////////////////////////////////////////
   while (!feof(f))
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
	  fgets(str2,100,f);
	
	///////////////////////////
	// if it's not a comment
	///////////////////////////
	else
	{
	  
	  //////////////////////////////////
	  // determine what is the input
	  /////////////////////////////////


  
	  

///////////////////////////////////////////////////////////////////////////
//			MOLECULE
///////////////////////////////////////////////////////////////////////////

	    ////////////////////////
	    // position of the atoms
	    ///////////////////////    
	     if( !strcmp(str,"position") || !strcmp(str,"POSITIONS")   || 
		!strcmp(str,"pos")  || !strcmp(str,"POS")  ||
	        !strcmp(str,"molecule") || !strcmp(str,"MOLECULE")	)
	    {  
		  fscanf(f,"%s",str_out);		  
		  strcpy(sys_xyz->pos,str_out);
	    }
	    
	    
	    ////////////////////////
	    // method to comoute H
	    ///////////////////////    
	    else if( !strcmp(str,"elec_struct") || !strcmp(str,"ELEC_STRUCT")  )
	    {  
		  fscanf(f,"%s",str_out);		  
		  strcpy(sys_xyz->elec_struct,str_out);
	    }
	    
	     ////////////////////////
	    // use S or not
	    ///////////////////////    
	    else if( !strcmp(str,"use_overlap") || !strcmp(str,"USE_OVERLAP")  )
	    {  
		  fscanf(f,"%s",str_out);		  
		  strcpy(sys_xyz->use_overlap,str_out);
	    }

	    /////////////////////////////////////////////////////
	    // orbitale electrode (USE ALL THE ORBITAL OR NOT)
	    /////////////////////////////////////////////////////    
	     else if( !strcmp(str,"orb_elec") || !strcmp(str,"ORB_ELEC")   || 
		!strcmp(str,"orb_cluster")  ||  !strcmp(str,"ORB_CLUSTER")  ||
	        !strcmp(str,"clust_orb") || !strcmp(str,"CLUST_ORB")	)
	    {  
		  fscanf(f,"%s",str_out);		  
		  strcpy(sys_xyz->orb_elec,str_out);
	    }
	    
	    
	    
	    ////////////////////////
	    // export  MO
	    ///////////////////////    
	     else if( !strcmp(str,"export_mo") || !strcmp(str,"EXPORT_MO")   )
	    {  
		  fscanf(f,"%s",str_out);		  
		  strcpy(sys_xyz->export_MO,str_out);
	    }
	   
    
      ////////////////////////
      // Number MO1 
      ///////////////////////
	     else if( !strcmp(str,"number_mo_1") )
       {
         fscanf(f,"%d",&INT);
         sys_xyz->nbr_mo_1=INT;
         sys_xyz->index_mo_1 = calloc(INT,sizeof(int));
       }
    
    
      ////////////////////////
      // index mo1
      ///////////////////////
	     else if( !strcmp(str,"index_mo_1") )
       {
         for(i=0;i<sys_xyz->nbr_mo_1;i++)
         {
           fscanf(f,"%d",&INT);
           sys_xyz->index_mo_1[i] = INT;
         }
       }

    
    
      ////////////////////////
      // Number MO2
      ///////////////////////
	     else if( !strcmp(str,"number_mo_2") )
       {
         fscanf(f,"%d",&INT);
         sys_xyz->nbr_mo_2=INT;
         sys_xyz->index_mo_2 = calloc(INT,sizeof(int));
       }
    
      ////////////////////////
      // index mo2
      ///////////////////////
	     else if( !strcmp(str,"index_mo_2") )
       {
         for(i=0;i<sys_xyz->nbr_mo_2;i++)
         {
           fscanf(f,"%d",&INT);
           sys_xyz->index_mo_2[i] = INT;
         }
       }
    
    
    ////////////////////////
    // Number MO3
    ///////////////////////
	     else if( !strcmp(str,"number_mo_3") )
       {
         fscanf(f,"%d",&INT);
         sys_xyz->nbr_mo_3=INT;
         sys_xyz->index_mo_3 = calloc(INT,sizeof(int));
       }
    
    ////////////////////////
    // index mo2
    ///////////////////////
	     else if( !strcmp(str,"index_mo_3") )
       {
         for(i=0;i<sys_xyz->nbr_mo_3;i++)
         {
           fscanf(f,"%d",&INT);
           sys_xyz->index_mo_3[i] = INT;
         }
       }
    
    
////////////////////////////////////////////////////////////////////////////////
// 				default case
////////////////////////////////////////////////////////////////////////////////
		   
	  else
	  {
		printf(" ERROR: CORRUPTED MAIN INPUT FILE: %s is not a valid entry \n",str);
		exit(1);
	  }
	    
	}
   }
   
   //close file
   fclose(f);
}