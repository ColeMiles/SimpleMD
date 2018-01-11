CXX=g++
CXXFLAGS=-std=c++11 -O2 -g

driver: driver.o simpleMD.o
	$(CXX) $(CXXFLAGS) -o driver driver.o simpleMD.o

driver.o: driver.cpp

simpleMD.o: simpleMD.cpp

clean:
	rm *.o
