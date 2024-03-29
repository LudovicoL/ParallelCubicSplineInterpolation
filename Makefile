.DEFAULT_GOAL := run
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
PRERUNFLAGS = -np 2
POSTRUNFLAGS = -fi "./input.txt" -s 0.1 -fo "./output.txt" -o n
HELPFLAGS = -help


install:
	sudo apt install build-essential
	sudo apt install mpich

compile:
	gcc ./gen_vectors.c -o ./gen
	gcc $(PREFLAGS) ./scsi.c -o ./scsi $(POSTFLAGS)
	mpicc $(PREFLAGS) $(FUNC) $(MAINC) -o $(MAINO) $(POSTFLAGS)
	
run: 
	mpirun $(PRERUNFLAGS) $(MAINO) $(POSTRUNFLAGS)

srun:
	./scsi $(POSTRUNFLAGS)

generate:
	./gen 100000 "./input.txt"

help:
	mpirun $(PRERUNFLAGS) $(MAINO) $(HELPFLAGS)
