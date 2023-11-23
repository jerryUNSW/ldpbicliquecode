CC=g++
CPPFLAGS=-I.
LDFLAGS=-g -fopenmp
DEPS = bigraph.h utility.h ldp-btf.h
OBJ = bigraph.o main.o utility.o ldp-btf.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)  

ldp-btf: $(OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) -lgomp 

clean:
	-rm -f ldp-btf *.o