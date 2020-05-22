CPP = g++

OPT = -c -O3 -Wall -Werror -m64 -std=c++11

OBJECTS_FTB = FTB.o

GSLLIB = -lgsl -lgslcblas -lm

#example
example.exe: example.o $(OBJECTS_FTB)
	$(CPP) -o example.exe example.o $(OBJECTS_FTB) $(GSLLIB)


FTB.o:	FBT.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ FBT.cpp
	@echo "..................done FBT."

example.o: example.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ example.cpp
	@echo "..................done Example."

all:
	make example.exe

clean:
	rm *.o
	rm *.exe
