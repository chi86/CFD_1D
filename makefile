# the compiler: gcc for C program, define as g++ for C++
CC  = gcc
CPP = g++
MAKEFLAGS = --silent

# special include directories
INCLUDE = -I.

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# laminar
# CFLAGS += -DLAMINAR

# if commented out --> singel precision
CFLAGS += -DUSE_DOUBLES=1

DEPS = main.o Cell.o Mesh.o Utilities.o CFD.o kEv2f.o

# temperatur dependent fluid properties
all: pipe_var
# constant material fluid properties
#all: pipe
#all: channel


verbose: CFLAGS += -DVERBOSE -g -Warray-bounds # -dM -E
verbose: all

debug: CFLAGS += -DDEBUG -g -Warray-bounds -fsanitize=address # -dM -E
debug: all

pipe: $(DEPS) Pipe.o 
	$(CC) $(CFLAGS) -o jasmine $+ -lm

# pipe_var: $(DEPS) Fluid_oil.o
# 	$(CC) $(CFLAGS) -o $@ $^ -lm

pipe_var: CFLAGS += -DVARMATPROPS
pipe_var: $(DEPS) Pipe.o Fluid_oil.o
	$(CC) $(CFLAGS) -o CFD_1D $^ -lm

channel: $(DEPS) Channel.o 
	$(CC) $(CFLAGS) -o CFD_1D $+ -lm

%.o: %.c
	$(CC) $(CFLAGS) $< -c


.PHONY: all

clean:
	$(RM) pipe *.o *~
