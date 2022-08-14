# Makefile for Jacobi_diagonalization

SHELL := cmd

CC := gcc
CLINKER := $(CC)

OPTS ?= -g
OPTLV ?= -O0

LAPACKROOT ?= C:\lapack-3.10.1
LIBPATH ?= -L $(LAPACKROOT)\lib
INCPATH ?= -I $(LAPACKROOT)\include
LIB ?= -l lapacke -l lapack -l cblas -l blas

.PHONY: all
all: Jacobi_diagonalization.exe

Jacobi_diagonalization.exe: Jacobi_diagonalization.obj compare.obj
	@echo Linking $@ against $^ and lapacke ...
	$(CLINKER) -o $@ $^ -static $(LIBPATH) $(LIB) -l gfortran -l quadmath $(OPTS)

Jacobi_diagonalization.obj: Jacobi_diagonalization.c
	@echo Compiling $@ from $< ...
	$(CC) -o $@ -c $< $(OPTLV) $(OPTS)

compare.obj: compare.c
	@echo Compiling $@ from $< ...
	$(CC) -o $@ -c $< $(OPTLV) $(INCPATH) $(OPTS)

.PHONY: veryclean
veryclean: clean
	-del /q Jacobi_diagonalization.exe 2> NUL

.PHONY: clean
clean:
	-del /q Jacobi_diagonalization.obj compare.obj 2> NUL

