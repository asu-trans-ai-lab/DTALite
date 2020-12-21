CC := g++
LIB := libstalite.so

all: $(LIB)

MathLibrary.o: MathLibrary.cpp
	$(CC) -c -fopenmp -fpic MathLibrary.cpp

$(LIB): MathLibrary.o
	$(CC) -shared -o $(LIB) MathLibrary.o

clean:
	rm -f *.o *.so
