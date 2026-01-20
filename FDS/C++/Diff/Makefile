CXX = g++
CXXFLAGS = -O2 -march=native
INCDIR = ../common

OBJ = main.o dif1d.o $(INCDIR)/bound.o
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

dif1d.o : dif1d.hpp

main.o : dif1d.o dif1d.hpp cnst.hpp $(INCDIR)
