CC=g++
CPPFLAGS=-I.
LDFLAGS=-g -fopenmp -lsqlite3
DEPS = bigraph.h utility.h biclique.h
OBJ = bigraph.o main.o utility.o biclique.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)  

biclique: $(OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) -lgomp 

testp3: test_p3_batch_with_ground_truth.cpp biclique.cpp bigraph.cpp utility.cpp
	$(CC) -std=c++17 -O3 -fopenmp -I. -o test_p3_batch_with_ground_truth test_p3_batch_with_ground_truth.cpp biclique.cpp bigraph.cpp utility.cpp -lsqlite3

clean:
	-rm -f biclique test_p3_batch_with_ground_truth *.o
