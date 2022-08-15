# Makefile for Jacobi_diagonalization

SHELL := cmd

FC := gfortran
FLINKER := $(FC)

OPTS ?= -s
OPTLV ?= -O2

LAPACKROOT ?= C:\lapack-3.10.1
LIBPATH ?= -L $(LAPACKROOT)\lib
LIB ?= -l lapack -l blas

.PHONY: all
all: Jacobi_diagonalization.exe

Jacobi_diagonalization.exe: Jacobi_diagonalization.obj compare.obj
	@echo Linking $@ against $^ and lapack ...
	$(FLINKER) -o $@ $^ -static $(LIBPATH) $(LIB) $(OPTS)

%.obj: %.f90
	@echo Compiling $@ from $< ...
	$(FC) -o $@ -c $< $(OPTLV) $(OPTS)

.PHONY: veryclean
veryclean: clean
	-del /q Jacobi_diagonalization.exe 2> NUL

.PHONY: clean
clean:
	-del /q Jacobi_diagonalization.obj compare.obj 2> NUL

