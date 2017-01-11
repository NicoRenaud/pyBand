#include "../defMacro.h"
#include "./header.h"




////////////////////////////////////////
//	find index from the reduced PARAM LIST
////////////////////////////////////////

int findIndex_reduced(char *name, char *PARAM)
{
  
  int index = 0;
	FILE *Fparameter;
	int i,j;
	char atom_name[3];

	// open the parameter file and read it
	Fparameter=fopen(PARAM,"r");
  
	if(Fparameter == NULL)
	{
	  //printf("Huckel parameters not found at:  /home/nico/Bureau/husky.1/SRC/huckel/parameters \n");
    printf("Huckel parameters not found at:  %s \n",PARAM);
	  exit(1);
	}


	for(i=1;i<=250;i++)
	{
		
		// read the atom name
		fscanf(Fparameter,"%s",atom_name);
		

		
		if (strcmp(name,atom_name) == 0)
		{
			index = i;
			break;
		}
		
		// read the rest
		fscanf(Fparameter,"%*i");
	
		for(j=0;j<4;j++)
		{
			fscanf(Fparameter,"%*i ");
		}
	
		for(j=0;j<4;j++)
			fscanf(Fparameter,"%*f %*f %*f %*f %*f %*f");
		fscanf(Fparameter,"\n");

  }
	fclose(Fparameter);

	// return the index 
  return(index);

}



////////////////////////////////////////
//	find index from the total PARAM LIST
////////////////////////////////////////

int findIndex(char *name)
{
  
  int index;
  
  if(!strcmp(name,"H"))	index = 1;
  else if(!strcmp(name,"He"))	index = 2;
  
  else if(!strcmp(name,"Li"))	index = 3;
  else if(!strcmp(name,"Be"))	index = 4;
  else if(!strcmp(name,"B"))	index = 5;
  else if(!strcmp(name,"C"))	index = 6;
  else if(!strcmp(name,"N"))	index = 7;
  else if(!strcmp(name,"O"))	index = 8;
  else if(!strcmp(name,"F"))	index = 9;
  else if(!strcmp(name,"Ne"))	index = 10;
  
  else if(!strcmp(name,"Na"))	index = 11;
  else if(!strcmp(name,"Mg"))	index = 12;
  else if(!strcmp(name,"Al"))	index = 13;
  else if(!strcmp(name,"Si"))	index = 14;
  else if(!strcmp(name,"P"))	index = 15;
  else if(!strcmp(name,"S"))	index = 16;
  else if(!strcmp(name,"Cl"))	index = 17;
  else if(!strcmp(name,"Ar"))	index = 18;
  
  else if(!strcmp(name,"K"))	index = 19;
  else if(!strcmp(name,"Ca"))	index = 20;
  else if(!strcmp(name,"Sc"))	index = 21;
  else if(!strcmp(name,"Ti"))	index = 22;
  else if(!strcmp(name,"V"))	index = 23;
  else if(!strcmp(name,"Cr"))	index = 24;
  else if(!strcmp(name,"Mn"))	index = 25;
  else if(!strcmp(name,"Fe"))	index = 26;
  else if(!strcmp(name,"Co"))	index = 27;
  else if(!strcmp(name,"Ni"))	index = 28;
  else if(!strcmp(name,"Cu"))	index = 29;
  else if(!strcmp(name,"Zn"))	index = 30;
  else if(!strcmp(name,"Ga"))	index = 31;
  else if(!strcmp(name,"Ge"))	index = 32;
  else if(!strcmp(name,"As"))	index = 33;
  else if(!strcmp(name,"Se"))	index = 34;
  else if(!strcmp(name,"Br"))	index = 35;
  else if(!strcmp(name,"Kr"))	index = 36;
  
  else if(!strcmp(name,"Rb"))	index = 37;
  else if(!strcmp(name,"Sr"))	index = 38;
  else if(!strcmp(name,"Y"))	index = 39;
  else if(!strcmp(name,"Zr"))	index = 40;
  else if(!strcmp(name,"Nb"))	index = 41;
  else if(!strcmp(name,"Mo"))	index = 42;
  else if(!strcmp(name,"Tc"))	index = 43;
  else if(!strcmp(name,"Ru"))	index = 44;
  else if(!strcmp(name,"Rh"))	index = 45;
  else if(!strcmp(name,"Pd"))	index = 46;
  else if(!strcmp(name,"Ag"))	index = 47;
  else if(!strcmp(name,"Cd"))	index = 48;
  else if(!strcmp(name,"In"))	index = 49;
  else if(!strcmp(name,"Sn"))	index = 50;
  else if(!strcmp(name,"Sb"))	index = 51;
  else if(!strcmp(name,"Te"))	index = 52;
  else if(!strcmp(name,"I"))	index = 53;
  else if(!strcmp(name,"Xe"))	index = 54;
  
  
  else if(!strcmp(name,"Cs"))	index = 55;
  else if(!strcmp(name,"Ba"))	index = 56;
  else if(!strcmp(name,"La"))	index = 57;
  else if(!strcmp(name,"Ce"))	index = 58;
  else if(!strcmp(name,"Pr"))	index = 59;
  else if(!strcmp(name,"Nd"))	index = 60;
  else if(!strcmp(name,"Pm"))	index = 61;
  else if(!strcmp(name,"Sm"))	index = 62;
  else if(!strcmp(name,"Eu"))	index = 63;
  else if(!strcmp(name,"Gd"))	index = 64;
  else if(!strcmp(name,"Tb"))	index = 65;
  else if(!strcmp(name,"Dy"))	index = 66;
  else if(!strcmp(name,"Ho"))	index = 67;
  else if(!strcmp(name,"Er"))	index = 68;
  else if(!strcmp(name,"Tm"))	index = 69;
  else if(!strcmp(name,"Yb"))	index = 70;
  
  
  else if(!strcmp(name,"Lu"))	index = 71;
  else if(!strcmp(name,"Hf"))	index = 72;
  else if(!strcmp(name,"Ta"))	index = 73;
  else if(!strcmp(name,"W"))	index = 74;
  else if(!strcmp(name,"Re"))	index = 75;
  else if(!strcmp(name,"Os"))	index = 76;
  else if(!strcmp(name,"Ir"))	index = 77;
  else if(!strcmp(name,"Pt"))	index = 78;
  else if(!strcmp(name,"Au"))	index = 79;
  else if(!strcmp(name,"Hg"))	index = 80;
  else if(!strcmp(name,"Ti"))	index = 81;
  else if(!strcmp(name,"Pb"))	index = 82;
  else if(!strcmp(name,"Bi"))	index = 83;
  else if(!strcmp(name,"Po"))	index = 84;
  else if(!strcmp(name,"At"))	index = 85;
  else if(!strcmp(name,"Rn"))	index = 86;
  
  else if(!strcmp(name,"Fr"))	index = 87;
  else if(!strcmp(name,"Ra"))	index = 88;
  else if(!strcmp(name,"Ac"))	index = 89;
  else if(!strcmp(name,"Th"))	index = 90;
  else if(!strcmp(name,"Pa"))	index = 91;
  else if(!strcmp(name,"U"))	index = 92;
  else if(!strcmp(name,"Np"))	index = 93;
  else if(!strcmp(name,"Pu"))	index = 94;
  else if(!strcmp(name,"Am"))	index = 95;
  else if(!strcmp(name,"Cm"))	index = 96;
  else if(!strcmp(name,"Bk"))	index = 97;
  else if(!strcmp(name,"Cf"))	index = 98;
  else if(!strcmp(name,"Es"))	index = 99;
  else if(!strcmp(name,"Fm"))	index = 100;
  else if(!strcmp(name,"Md"))	index = 101;
  else if(!strcmp(name,"No"))	index = 102;
  else if(!strcmp(name,"Lr"))	index = 103;
  
  else
  {
    //       printf("Element %s not recognized to find its index, please check XYZ__ file\n",name);
    index = 0;
    //       exit(0);
  }
  return(index);
  
}











