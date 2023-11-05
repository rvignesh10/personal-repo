CXX=mpic++
opt= -O3
CCFLAGS= -fPIC -O3 -I$(APlusPlus)/include -I$(PPlusPlus)/include
AppLibraries=-Wl,-rpath,$(PPlusPlus)/lib -L$(PPlusPlus)/lib -lPpp -lPpp_static
LIBS= $(AppLibraries)

heat2d:
	$(CXX) $(CCFLAGS) -c heat2d.C
	$(CXX) $(CCFLAGS) -o heat2d.bin heat2d.o
clean:
	rm *.o
	rm *.bin
	rm *.C~
	rm *.sh~
