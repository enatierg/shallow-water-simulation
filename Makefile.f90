#############################################################
# Makefile for Shallow Water Model code
#############################################################

# Compiler and flags
F95 = mpif90
CFLAGS = -g -Wall
LFLAGS = -g
# Link to BLAS library.
LTBS = -L${BLASDIR} ${BLASLIB} -lpthread

# Object files (except header compiled separately)
OBJS = matmult.o \
       eta_peak.o \
       create_matrices.o \
       leapfrog.o \
       save_field.o \
       average.o \
       sw_main.o \
       sparsegather.o \
       test_leapfrog.o

# Main executable
MAIN = sw_main

# Main target: Create executable
all: $(MAIN)

# Compile header first
header.mod: header.f90
	$(F95) -c $(CFLAGS) -o header.o header.f90

# Compile individual files
%.o: %.f90 header.mod
	$(F95) -c $(CFLAGS) -o $@ $<

# Link files to create main executable
$(MAIN): $(OBJS)
	$(F95) -o $(MAIN) $(LFLAGS) header.o $(OBJS) $(LTBS)

# Tidy up
clean:
	rm -f *.o *.mod core.* $(MAIN)
