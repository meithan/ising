# ==============================================================================
# Ising Model Makefile
#
# Available build targets:
#  'ising' (default): performs ising run(s) at a fixed temperature
#  'clean': removes all object files and the compiled binary
# ==============================================================================

# Compiler
#COMPILER= ifort
COMPILER= g++

# Additional user compilation flags
#USER_FLAGS = -g -Wall -pedantic
USER_FLAGS= -O3

# ==============================================================================

CFLAGS= $(USER_FLAGS)
PROGRAMS= ising

# ==============================================================================
# BUILD TARGETS

default : ising

ising : IsingModel.o ising.o
	$(COMPILER) $(CFLAGS) IsingModel.o ising.o -o ising

.PHONY: clean
clean :
	rm -f *.o $(PROGRAMS)

# ==============================================================================
# OBJECT BUILD RULES

IsingModel.o : IsingModel.cpp IsingModel.h utils.h
	$(COMPILER) $(CFLAGS) -c IsingModel.cpp

ising.o : IsingModel.cpp IsingModel.h utils.h ising.cpp
	$(COMPILER) $(CFLAGS) -c ising.cpp
