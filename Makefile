# Makefile for TwoPunctures

BASE=$(shell /bin/pwd)
SRCD=$(BASE)/src
INCD=$(BASE)/include
OBJD=$(BASE)/obj
LIBD=$(BASE)/lib


NAME=TwoPuncturesRun
LIBNAME=TwoPunctures
EXE=$(NAME).x
SRC=$(wildcard $(SRCD)/*.c)
OBJ=$(patsubst $(SRCD)/%.c,$(OBJD)/%.o,$(SRC))
INC=$(wildcard $(INCD)/*.c)
INC_PARAMS=$(foreach d, $(INCD), -I$d)

OBJ2=$(wildcard $(OBJD)/*.o)
LIB=$(LIBD)/lib$(LIBNAME).so


# mandatory flags

CC = gcc
LD = ld
AR = ar

CFLAGS = -std=c11 -fPIC -pedantic $(shell gsl-config --cflags)
LFLAGS = $(shell gsl-config --libs)

CFLAGS += -O3
CFLAGS += -Wall -Wwrite-strings -Wshadow
CFLAGS += -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -fshort-enums -fno-common

# old flag information [Leave alone]
##CFLAGS + =`gsl-config --cflags`# GSL
##CFLAGS += -fopenmp # activate OMP (opt)
##LFLAGS=`gsl-config --libs`


all: $(EXE) $(LIB)
	@echo "All done"

$(EXE): $(OBJ)
	@echo "Building $@ ..."
	@echo

	$(CC) $(CFLAGS) $(OBJ) $(LFLAGS) -o $@
        
#	# this doesn't need to be a shared object
#	@rm $(OBJD)/TwoPuncturesRun.o > /dev/null 2>&1

$(OBJD)/%.o: $(SRCD)/%.c
	@echo "Building objects ..."
	@echo

	@mkdir -p $(dir $@)
#	@echo "Compiling $< ... -> $@"
	$(CC) $(CFLAGS) $(INC_PARAMS) -c $< -o $@

$(LIB): $(OBJ)
	@echo "Making libraries... "
	@echo

	@mkdir -p $(LIBD)
	$(CC) $(CFLAGS) $(LFLAGS) -shared $(OBJ) -o $(LIBD)/lib$(LIBNAME).so

	$(AR) rcs $(LIBD)/libTwoPunctures.a $(OBJ)
#	$(AR) rcs $(LIBD)/libTwoPunctures_static.a $@ $^

clean:
	@echo "Cleaning ..."
	@rm -rf $(OBJD)
	@rm -rf $(LIBD)
	@rm -rf $(EXE)
	@echo "... done"

.PHONY: all clean
