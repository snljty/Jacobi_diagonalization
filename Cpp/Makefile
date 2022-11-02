# Makefile for Jacobi_diagonalization

SHELL = cmd
CXX := g++
CXXLINKER := $(CXX)

DEBUG_LEVEL := 

LAPACKROOT := C:\lapack-3.10.1
INCPATH := -I $(LAPACKROOT)\include
LIBPATH := -L $(LAPACKROOT)\lib
LIB := -l lapacke -l lapack -l cblas -l blas

ifneq ("$(DEBUG_LEVEL)","")
	OPTS := -g
	OPTLV := -O0
else
	OPTS := -D NDEBUG -s
	OPTLV := -O2
endif

.PHONY: all
all: test.exe

test.exe: test.obj Jacobi_diagonalization.obj
	@echo Linking $@ against $^ ...
	$(CXXLINKER) -o $@ $^ -static $(OPTS) $(LIBPATH) $(LIB) -l gfortran -l quadmath

Jacobi_diagonalization.obj: Jacobi_diagonalization.cpp Jacobi_diagonalization.hpp
	@echo Compiling $@ from $< ...
	$(CXX) -o $@ -c $< $(OPTLV) $(OPTS) $(INCPATH)

test.obj: test.cpp Jacobi_diagonalization.hpp
	@echo Compiling $@ from $< ...
	$(CXX) -o $@ -c $< $(OPTLV) $(OPTS)

.PHONY: clean
clean:
	-del /q test.obj Jacobi_diagonalization.obj 2> NUL

.PHONY: 
veryclean: clean
	-del /q test.exe 2> NUL

