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
sequence_t **count_and_read_seqs( char *filename );
static void substring_indices( char *src, char *dest, const int start, const int end );
static inline int num_substrings( const int str_len, const int window_size );

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

    sequence_t **refseqs     = NULL;
    sequence_t **design_seqs = NULL;

    refseqs     = count_and_read_seqs( ref_file_name );
    design_seqs = count_and_read_seqs( design_file_name );

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

sequence_t **count_and_read_seqs( char *filename )
{

    FILE *open_file = NULL;
    sequence_t **local_seqs = NULL;
    int num_seqs = 0;

    char local_filename[ MAX_STRING_SIZE ];

    strcpy( local_filename, filename );

    open_file = fopen( local_filename, "r" );

    num_seqs = count_seqs_in_file( open_file );

    local_seqs = malloc( sizeof( sequence_t *) * num_seqs );

    read_sequences( open_file, local_seqs );

    fclose( open_file );

    return local_seqs;

}

static void substring_indices( char *src, char *dest, const int start, const int end )
{
    int index       = 0;
    int inner_index = 0;

    for( index = start; index < end; index++ )
        {
            dest[ inner_index ] = src[ index ];

            dest[ ++inner_index ] = '\0';
        }
}

static inline int num_substrings( const int str_len, const int window_size )
{
    return str_len - window_size + 1;
}
