# Pre-options to compile the files:
PREFLAGS = -O3
# Post-options to compile the files:
POSTFLAGS = -lm

# External library
FUNC = ./MyLibrary.c

# Source codes and executables:
MAINC = ./pcsi.c
MAINO = ./pcsi

# Flag to the run:
PRERUNFLAGS = -np 3
POSTRUNFLAGS = -fi "./input.txt" -s 0.1 -fo "./output.txt"


install:
	sudo apt install mpich

compile:
	gcc ./gen_vectors.c -o ./gen
	mpicc $(PREFLAGS) $(FUNC) $(MAINC) -o $(MAINO) $(POSTFLAGS)
	
run: 
	mpirun $(PRERUNFLAGS) $(MAINO) $(POSTRUNFLAGS)

generate:
	./gen 10 "./input.txt"
