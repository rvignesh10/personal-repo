CXX=mpic++
opt= -O3
CCFLAGS= -fPIC -O3 

heat2d:
	$(CXX) $(CCFLAGS) -c heat2d.C
	$(CXX) $(CCFLAGS) -o heat2d.bin heat2d.o
clean:
	rm *.o
	rm *.bin
