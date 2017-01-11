#include "../defMacro.h"

////////////////////////////////////
// Print the header of the program //
////////////////////////////////////



void print_header_file_xyz(char *file_name,SYS_XYZ sys) {
 
  FILE *f;

  f = fopen(file_name,"w");
  if (!f) 
  {
      printf("can't create %s \n",file_name);
      exit(1);
  }
  fprintf(f,"\n====  ===============================  ====\n");    
  fprintf(f,"====                                   ====\n");   
  fprintf(f,"====  	    Huckel :)	       ====\n"); 
  fprintf(f,"====                                   ====\n");   
  fprintf(f,"====  ===============================  ====\n\n");
  
  
  fprintf(f," == Coordinates of the molecule    \t: %s\n",sys.pos);  
//   fprintf(f," == Export transmission            \t: %s\n",sys.export_TE);
//   fprintf(f," == Export molecular orbitals      \t: %s\n",sys.export_MO);
//   fprintf(f," == Export Current               \t: %s\n",sys.export_current);
  fclose(f);
}





