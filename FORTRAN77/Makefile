COMP = gfortran
FFlAGS = -O
SOURCE = FBT.f example.f

example:
	$(COMP) $(FFLAGS) -c $(SOURCE)

	$(COMP) *.o -o example.out

clean:
	rm *.o
	rm example.out
