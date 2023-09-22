# 
# THIS MAKEFILE USES THE FOLLOWING ENVIRONMENTAL VARIABLE:
#      APlusPlus
# The APlusPlus env variable should normally be set when you login.
#
# If not you can set it using:
#
# cg6: use:
# 	setenv APlusPlus /home/henshw/software/AppPpp-0.8.3-gcc7.5.0/A++/install
#
# cgpi: use:
# 	setenv APlusPlus /home/henshw/software/AppPpp-0.8.3-gcc4.8.5/A++/install
#
# If you type:
# 		ls $APlusPlus/include
# you should see A++.h as well as other include files.

# First target is made by default when using "make", traditionally named "all"
all = heat1dImp

# Set names of compilers once in case we need to change them
CC = gcc
CXX = g++

CCFLAGS= -fPIC -O3 -I$(APlusPlus)/include -I../../include

# List of libraries for A++
AppLibrariers = -Wl, -rpath, $(APlusPlus)/lib -L$(APlusPlus)/lib -lApp -lApp_static

# List of all libraries
LIBS = $(AppLibrariers)

# implicit (generic) rule to compile .C files
# $@ = file name of the target
# $< = name of the first prerequisite
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

# 1D heat equation, implicit time-stepping, A++ arrays
heat1dImpFiles = heat1dImp.o tridiagonal.o
heat1dImp: $(heat1dImpFiles)
	$(CXX) $(CCFLAGS) -o $@ $(heat1dImpFiles) $(LIBS)

clean: ; rm *.o