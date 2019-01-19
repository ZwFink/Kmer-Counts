#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <stdbool.h>

#include "hash_table.h"
#include "array_list.h"
#include "protein_oligo_library.h"

const int NUM_ARGS = 5;
const int MAX_STRING_SIZE = 512;

static inline bool tolerable_match( char *a, char *b, int size, int num_mismatches );

int main( int argc, char **argv )
{

    char design_file_name[ MAX_STRING_SIZE ];
    char ref_file_name[ MAX_STRING_SIZE ];
    char outfile_name[ MAX_STRING_SIZE ];

    int num_threads = 0;

    if( argc != NUM_ARGS )
        {
            printf( "USAGE: get_kmer_counts design_file_name "
                    "ref_file_name outfile_name num_threads\n"
                  );
            return EXIT_FAILURE;
        }

    strcpy( design_file_name, argv[ 1 ] );
    strcpy( ref_file_name,    argv[ 2 ] );
    strcpy( outfile_name,     argv[ 3 ] );

    num_threads = atoi( argv[ 4 ] );

    omp_set_num_threads( num_threads );

    return EXIT_SUCCESS;
}

static inline bool tolerable_match( char *a, char *b, int size, int num_mismatches )
{
    int index = 0;
    int mismatches = 0;

    if( a[ 0 ] == b[ 0 ]
        || a[ size ] == b[ size ]
      )
        {
            for( index = 0; index < size; index++ )
                {
                    if( a[ index ] - b[ index ] )
                        {
                            mismatches++;
                            if( mismatches > num_mismatches )
                                {
                                    return false;
                                }
                        }
                }
        }
    return true;
}
