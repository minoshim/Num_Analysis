CXX = g++
CXXFLAGS = -O2 -march=native

OBJ = bound.o

.PHONY : all clean

all : $(OBJ)

%.o : %.cpp
	$(CXX) -c $< $(CXXFLAGS)

clean :
	rm -f $(OBJ)

bound.o : bound.hpp
