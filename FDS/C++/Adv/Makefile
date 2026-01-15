CC = g++
CFLAGS = -O2 -march=native

OBJ = main.o adv1d.o bound.o

.PHONY : all
.PHONY : clean

all : $(OBJ)
	$(CC) -o a.out $(OBJ) $(CFLAGS)

%.o : %.cpp
	$(CC) -c $< $(CFLAGS)

clean :
	rm *.o *.out -f
cdata :
	rm *.dat -f

adv1d.o : adv1d.hpp
bound.o : bound.hpp

main.o : adv1d.o bound.o adv1d.hpp bound.hpp cnst.hpp routine.hpp
