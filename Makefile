# Makefile for TwoPunctures

BASE=$(shell /bin/pwd)
SRCD=$(BASE)/src
OBJD=$(BASE)/obj
LIBD=$(BASE)/lib

NAME=TwoPuncturesRun
EXE=$(NAME).x
SRC=$(wildcard $(SRCD)/*.c)
OBJ=$(patsubst $(SRCD)/%.c,$(OBJD)/%.o,$(SRC))
LIB=$(LIBD)/lib$(NAME).a

CC=gcc
CFLAGS=-g -std=c99 `gsl-config --cflags`
#CFLAGS=-O3 -std=c99 `gsl-config --cflags`
LFLAGS=`gsl-config --libs`
LDLFLAGS=-lgsl -lgslcblas -lm

all: $(EXE) # $(LIB)
	@echo "All done"

$(EXE): $(OBJ)
	@echo "Building $@ ..."
	@$(CC) $(CFLAGS) $(LFLAGS) $(OBJ) -o $@ $(LDLFLAGS)
        # this doesn't need to be a shared object
	rm $(OBJD)/TwoPuncturesRun.o > /dev/null 2>&1

$(OBJD)/%.o: $(SRCD)/%.c
	@mkdir -p $(dir $@)
	@echo "Compiling $< ... -> $@"
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	@echo "Cleaning ..."
	@rm -vf $(OBJ)
	@rm -vf $(LIB)
	@echo "... done"

new: clean $(EXE)
	@echo "All done"

clearup: clean
	@echo "Removing dirs and exe ..."
	@rm -rvf $(OBJ)
	@rm -rvf $(LIB)
	@rm -rvf $(EXE)
	@echo "All done"

.PHONY: all clean new #help

# DO NOT DELETE
