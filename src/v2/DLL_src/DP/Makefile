CC := g++
LIB := libstalite.so

all: $(LIB)

main_api.o : main_api.cpp
	$(CC) -c -fopenmp -fpic main_api.cpp

$(LIB): main_api.o
	$(CC) -shared -o $(LIB) main_api.o

clean:
	rm -f *.o *.so
