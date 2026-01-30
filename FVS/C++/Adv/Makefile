CXX = g++
CXXFLAGS = -O2 -march=native
INCDIR = ../common

OBJ = main.o adv1d.o $(INCDIR)/bound.o
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

adv1d.o : adv1d.hpp

main.o : adv1d.o adv1d.hpp cnst.hpp $(INCDIR)
