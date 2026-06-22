INCDIR = ../common
include $(INCDIR)/Makefile.inc

OBJ = main.o bgs1d.o $(INCDIR)/bound.o
TARGET = a.out
COMMON_HEADERS = $(INCDIR)/routine.hpp $(INCDIR)/bound.hpp

.PHONY : all clean cdata

all : $(TARGET)

$(TARGET) : $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

%.o : %.cpp
	$(CXX) -c $< $(CXXFLAGS) -I$(INCDIR)

$(INCDIR)/bound.o: $(INCDIR)/bound.cpp $(INCDIR)/bound.hpp $(INCDIR)/Makefile $(INCDIR)/Makefile.inc
	$(MAKE) -C $(INCDIR)

clean :
	rm -f *.o $(TARGET)
cdata :
	rm -f *.dat

bgs1d.o : bgs1d.hpp

main.o : bgs1d.o bgs1d.hpp cnst.hpp $(COMMON_HEADERS)
