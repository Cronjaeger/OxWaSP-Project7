# Makefile for building the C language shared libraries for paralellComparisonpkg-package.
C = gcc
OPTS = -c -fPIC
LOADER = gcc

OBJECTS = square.o

 square.so: $(OBJECTS)
  R CMD SHLIB -o square.so $(OBJECTS)

.c.o: ; $(C) $(OPTS) -c $<

clean:
	-rm *.o *.so