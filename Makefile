CFLAGS= -O0 -Wall -Wextra -std=c99 -pedantic -lpthread

get_kmer_counts: get_kmer_counts.o protein_oligo_library.o dynamic_string.o hash_table.o array_list.o set.o 
	gcc $(CFLAGS) get_kmer_counts.o protein_oligo_library.o dynamic_string.o hash_table.o array_list.o set.o -o get_kmer_counts 
get_kmer_counts.o: get_kmer_counts.c protein_oligo_library.h hash_table.h array_list.h set.h

protein_oligo_library.o: protein_oligo_library.c protein_oligo_library.h hash_table.h array_list.h set.h

dynamic_string.o: dynamic_string.c dynamic_string.h

array_list.o: array_list.c array_list.h 

hash_table.o: hash_table.c hash_table.h

array_list: array_list_main.o array_list.o
array_list_main.o: array_list_main.c array_list.h
array_list.o: array_list.c array_list.h

set.o: set.c set.h 


.PHONY: debug clean optimized profile
debug: CFLAGS+= -g -O0 
debug: clean
debug: get_kmer_counts

optimized: CFLAGS += -O3  -ffast-math -fopenmp
optimized: clean
optimized: get_kmer_counts

profile: CFLAGS += -pg -g
profile: clean
profile: get_kmer_counts


clean:
	rm -rf *.o *.gch get_kmer_counts

