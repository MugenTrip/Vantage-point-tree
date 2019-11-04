# define the shell to bash
SHELL := /bin/bash

IDIR =./inc

ODIR=./src

SDIR=./src

_OBJ = vptree_sequential.o vptree_pthreads.o vptree_openmp.o vptree_cilk.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_SRC = vptree_sequential.c vptree_pthreads.c vptree_openmp.c vptree_cilk.c
SRC = $(patsubst %,$(SDIR)/%,$(_SRC))

_DEPS = quickselect.h vptree.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_FINAL = vptree_sequential.a vptree_pthreads.a vptree_openmp.a vptree_cilk.a
FINAL = $(patsubst %,$(IDIR)/%,$(_FINAL))

# define the C/C++ compiler to use,default here is clang
CC = gcc-5

CFLAGS = -O3 -fopenmp -fcilkplus -Wall -lpthread -I$(IDIR)

all: $(FINAL)

$(IDIR)/%.a: $(OBJ)
	ar rcs inc/vptree_sequential.a src/vptree_sequential.o
	ar rcs inc/vptree_pthreads.a  src/vptree_pthreads.o
	ar rcs inc/vptree_openmp.a  src/vptree_openmp.o
	ar rcs inc/vptree_cilk.a  src/vptree_cilk.o

$(ODIR)/%.c.o: $(SRC) $(DEPS)
	$(CC) -c $@ -o $(CFLAGS)




