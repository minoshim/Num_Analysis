INCDIR = ../common
include $(INCDIR)/Makefile.inc

OBJ = main.o dif1d.o $(INCDIR)/bound.o
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

dif1d.o : dif1d.hpp

main.o : dif1d.o dif1d.hpp cnst.hpp $(COMMON_HEADERS)
