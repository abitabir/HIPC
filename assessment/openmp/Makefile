SHELL:=/bin/bash
CC=gcc
MM=module
CFLAGS=-O3 
MATHFLAG=-lm
PAPIFLAG=-lpapi
PGFLAG=-pg
OPENMPFLAG=-fopenmp

OBJDIR = obj

_OBJ = main.o args.o vtk.o data.o setup.o boundary.o vortex.o serial.o data_serial.o setup_serial.o boundary_serial.o vortex_serial.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

.PHONY: directories

all: directories vortex

start: libs directories vortex

obj/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

libs:
	$(MM) load perf/PAPI/5.6.0-GCCcore-7.3.0
	$(MM) load compiler/GCC/11.2.0

vortex: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(MATHFLAG) $(PAPIFLAG) $(OPENMPFLAG) $(PGFLAG)

clean:
	rm -Rf $(OBJDIR)
	rm -f vortex

directories: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

