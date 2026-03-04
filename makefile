CC=g++
CPPFLAGS=-I.
LDFLAGS=-g -fopenmp
DEPS = bigraph.h utility.h biclique.h
OBJ = bigraph.o main.o utility.o biclique.o sqlite3_c.o

# Compile sqlite3.c as C (not C++)
sqlite3_c.o: sqlite3.c
	gcc -c -O2 -o sqlite3_c.o sqlite3.c -DSQLITE_THREADSAFE=1

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)  

biclique: $(OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) -lgomp -ldl

# Exact DDS solver
exact_dds: exact_dds.cpp biclique.cpp bigraph.cpp utility.cpp sqlite3_c.o
	$(CC) -std=c++17 -O3 -fopenmp -I. -o exact_dds exact_dds.cpp biclique.cpp bigraph.cpp utility.cpp sqlite3_c.o -ldl -lpthread

# ELDP DDS algorithm
biclique_dds: main_dds.cpp biclique_dds_eldp.cpp biclique_dds_eldp_v12.cpp biclique_dds_eldp_v13.cpp biclique_dds_eldp_v14.cpp biclique_dds_eldp_v15.cpp biclique_dds_eldp_v16.cpp biclique_dds_eldp_v17.cpp biclique_dds_naive.cpp biclique_dds_naive_eldp_v2.cpp biclique_dds_peeling_eldp.cpp biclique_dds_peeling_expand_shrink_eldp.cpp biclique.cpp bigraph.cpp utility.cpp sqlite3_c.o
	$(CC) -std=c++17 -O3 -fopenmp -I. -o biclique_dds main_dds.cpp biclique_dds_eldp.cpp biclique_dds_eldp_v12.cpp biclique_dds_eldp_v13.cpp biclique_dds_eldp_v14.cpp biclique_dds_eldp_v15.cpp biclique_dds_eldp_v16.cpp biclique_dds_eldp_v17.cpp biclique_dds_naive.cpp biclique_dds_naive_eldp_v2.cpp biclique_dds_peeling_eldp.cpp biclique_dds_peeling_expand_shrink_eldp.cpp biclique.cpp bigraph.cpp utility.cpp sqlite3_c.o -ldl -lpthread

testp3: test_p3_batch_with_ground_truth.cpp biclique.cpp bigraph.cpp utility.cpp sqlite3_c.o
	$(CC) -std=c++17 -O3 -fopenmp -I. -o test_p3_batch_with_ground_truth test_p3_batch_with_ground_truth.cpp biclique.cpp bigraph.cpp utility.cpp sqlite3_c.o -ldl -lpthread

variance_analysis: variance_analysis.cpp bigraph.cpp utility.cpp sqlite3_c.o
	$(CC) -std=c++17 -O3 -fopenmp -I. -o variance_analysis variance_analysis.cpp bigraph.cpp utility.cpp sqlite3_c.o -ldl -lpthread

verify_unbiasedness: verify_unbiasedness.cpp biclique_dds_eldp.cpp biclique.cpp bigraph.cpp utility.cpp sqlite3_c.o
	$(CC) -std=c++17 -O3 -fopenmp -I. -o verify_unbiasedness verify_unbiasedness.cpp biclique_dds_eldp.cpp biclique.cpp bigraph.cpp utility.cpp sqlite3_c.o -ldl -lpthread

clean:
	-rm -f biclique exact_dds biclique_dds test_p3_batch_with_ground_truth variance_analysis *.o
