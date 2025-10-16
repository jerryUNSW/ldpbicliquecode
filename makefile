CC=g++
CPPFLAGS=-I.
LDFLAGS=-g -fopenmp -lsqlite3
DEPS = bigraph.h utility.h biclique.h
OBJ = bigraph.o main.o utility.o biclique.o
P3_OBJ = bigraph.o test_p3_batch_with_ground_truth.o utility.o biclique.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)  

biclique: $(OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) -lgomp 

test_p3_batch: $(P3_OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) -lgomp 

clean:
	-rm -f biclique test_p3_batch *.o
