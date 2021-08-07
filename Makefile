SHELL:='/bin/bash'
CXX=g++
IDIR=include/mcpm/
IDIRS=-Iinclude/
CFLAGS=-std=c++17 -Wall -O2 $(IDIRS) $(LDIR) -fopenmp
ODIR=obj/
SRCDIR=src/

LIB_PATH=lib/

_DEPS=graph.h binary_heap.h matching.h globals.h
DEPS=$(patsubst %,$(IDIR)%,$(_DEPS))

_OBJ_ALL=graph.o binary_heap.o matching.o
OBJ_ALL=$(patsubst %,$(ODIR)%,$(_OBJ_ALL))

all: $(LIB_PATH)libmcpm.a 

$(ODIR)%.o: $(SRCDIR)%.cc $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

$(LIB_PATH)libmcpm.a: $(OBJ_ALL)
	ar rcs $@ $^

bin/example: example.cc
	$(CXX) example.cc -Llib/ -lmcpm -Iinclude/ -o bin/example
	
.PHONY: clean

clean:
	rm $(ODIR)*.o $(LIB_PATH)*
