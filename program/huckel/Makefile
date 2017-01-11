# ==== ====================================== ====
# ==== 	        HUCKEL  Makefile 	      ====
# ==== ====================================== ====
# Automatically generated the  21-07-2014 at 11h55


CC   = gcc

OBJ = OBJ/algebra_cplx.o OBJ/algebra.o OBJ/extract_hamiltonian.o OBJ/huckel.o OBJ/read_input_hkl.o OBJ/read_input_hkl_general.o OBJ/main.o OBJ/create_system.o OBJ/read_input_file.o OBJ/export_gnuplot.o OBJ/print_file.o OBJ/print_screen.o 

LIBS = -lm -lgslcblas -llapack -lgsl
BIN  =  huckel
CFLAGS = -Wall
RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all: all-before $(BIN) all-after

clean: clean-custom
		${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $(BIN) $(LIBS)

OBJ/algebra_cplx.o : ./SRC/algebra_cplx/algebra_cplx.c
	$(CC) -c ./SRC/algebra_cplx/algebra_cplx.c -o OBJ/algebra_cplx.o $(CFLAGS)

OBJ/algebra.o : ./SRC/algebra_real/algebra.c
	$(CC) -c ./SRC/algebra_real/algebra.c -o OBJ/algebra.o $(CFLAGS)

OBJ/extract_hamiltonian.o : ./SRC/algebra_real/extract_hamiltonian.c
	$(CC) -c ./SRC/algebra_real/extract_hamiltonian.c -o OBJ/extract_hamiltonian.o $(CFLAGS)

OBJ/huckel.o : ./SRC/huckel/huckel.c
	$(CC) -c ./SRC/huckel/huckel.c -o OBJ/huckel.o $(CFLAGS)

OBJ/read_input_hkl.o : ./SRC/huckel/read_input_hkl.c
	$(CC) -c ./SRC/huckel/read_input_hkl.c -o OBJ/read_input_hkl.o $(CFLAGS)

OBJ/read_input_hkl_general.o : ./SRC/huckel/read_input_hkl_general.c
	$(CC) -c ./SRC/huckel/read_input_hkl_general.c -o OBJ/read_input_hkl_general.o $(CFLAGS)

OBJ/main.o : ./SRC/main.c
	$(CC) -c ./SRC/main.c -o OBJ/main.o $(CFLAGS)

OBJ/create_system.o : ./SRC/manag_input/create_system.c
	$(CC) -c ./SRC/manag_input/create_system.c -o OBJ/create_system.o $(CFLAGS)

OBJ/read_input_file.o : ./SRC/manag_input/read_input_file.c
	$(CC) -c ./SRC/manag_input/read_input_file.c -o OBJ/read_input_file.o $(CFLAGS)

OBJ/export_gnuplot.o : ./SRC/manag_output/export_gnuplot.c
	$(CC) -c ./SRC/manag_output/export_gnuplot.c -o OBJ/export_gnuplot.o $(CFLAGS)

OBJ/print_file.o : ./SRC/manag_output/print_file.c
	$(CC) -c ./SRC/manag_output/print_file.c -o OBJ/print_file.o $(CFLAGS)

OBJ/print_screen.o : ./SRC/manag_output/print_screen.c
	$(CC) -c ./SRC/manag_output/print_screen.c -o OBJ/print_screen.o $(CFLAGS)

