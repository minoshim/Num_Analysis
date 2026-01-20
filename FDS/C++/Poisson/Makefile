CXX = g++
CXXFLAGS = -O2 -march=native
INCDIR = ../common

OBJ = main.o poi2d.o $(INCDIR)/bound.o
TARGET = a.out

.PHONY : all clean cdata

all : $(TARGET)

$(TARGET) : $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

%.o : %.cpp
	$(CXX) -c $< $(CXXFLAGS) -I$(INCDIR)

$(INCDIR)/bound.o:
	make -C $(INCDIR)

clean :
	rm -f *.o $(TARGET)
cdata :
	rm -f *.dat

poi2d.o : poi2d.hpp

main.o : poi2d.o poi2d.hpp cnst.hpp routine.hpp $(INCDIR)
