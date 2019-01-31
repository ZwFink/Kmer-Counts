#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <stdbool.h>

#include "hash_table.h"
#include "array_list.h"
#include "protein_oligo_library.h"

#ifndef _OPENMP
    #define omp_get_wtime() 0
#endif

const int NUM_ARGS         = 5;
const int WINDOW_SIZE      = 9;
const int NUM_MISMATCHES   = 1;
const int MAX_STRING_SIZE  = 512;
const int LARGE_TABLE_SIZE = 4000000;

typedef struct kmer
{
    char *seq;
    unsigned int kmer_start;
    unsigned int kmer_end;
    unsigned int kmer_score;
} kmer_t;

static inline bool tolerable_match( char *a, char *b, int size, int num_mismatches );
static inline void add_valid_kmers( kmer_t **kmers, const unsigned int num_subsets, hash_table_t *table );
static inline void copy_kmer( kmer_t *dest, kmer_t *src );
sequence_t **count_and_read_seqs( char *filename );
static void substring_indices( char *src, char *dest, const int start, const int end );
static inline int num_substrings( const int str_len, const int window_size );
static void subset_lists_local( kmer_t **dest_arr, char *seq,
                                int sequence_len, const int window_size );
static void subset_lists_ht( hash_table_t *dest, char *seq,
                              int sequence_len, const int window_size );

hash_table_t *seqs_to_kmer_table( sequence_t **seqs, const int num_seqs );
void get_kmer_totals( hash_table_t *target_kmers, sequence_t **designed_oligos, int num_oligos, int num_mismatches );
void get_mismatch_counts( hash_table_t *table, HT_Entry **items, char *kmer,
                          unsigned int num_items, int num_mismatches
                        );
void write_outputs( char *out_file, hash_table_t *table );
void clear_table( hash_table_t *table );
void kmer_init( kmer_t *kmer, char *seq, unsigned int start, unsigned int end, unsigned int score );

int main( int argc, char **argv )
{

    char design_file_name[ MAX_STRING_SIZE ];
    char ref_file_name[ MAX_STRING_SIZE ];
    char outfile_name[ MAX_STRING_SIZE ];

    int num_seqs_ref    = 0;
    int num_seqs_design = 0;

    FILE *open_file = NULL;

    int num_threads = 0;

    hash_table_t *target_seqs = NULL;

    double start_time = 0;
    double end_time   = 0;

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

    open_file = fopen( ref_file_name, "r" );
    num_seqs_ref = count_seqs_in_file( open_file );
    fclose( open_file );

    open_file = fopen( design_file_name, "r" );
    num_seqs_design = count_seqs_in_file( open_file );
    fclose( open_file );

    num_threads = atoi( argv[ 4 ] );

    #ifdef _OPENMP
    omp_set_num_threads( num_threads );
    #endif

    sequence_t **refseqs     = NULL;
    sequence_t **design_seqs = NULL;

    start_time = omp_get_wtime();
    
    refseqs     = count_and_read_seqs( ref_file_name );
    design_seqs = count_and_read_seqs( design_file_name );

    target_seqs = seqs_to_kmer_table( refseqs, num_seqs_ref );

    get_kmer_totals( target_seqs, design_seqs,
                     num_seqs_design, NUM_MISMATCHES
                   );

    write_outputs( outfile_name, target_seqs );

    end_time = omp_get_wtime();

    printf( "Finished in %f seconds\n", end_time - start_time );

    clear_table( target_seqs );

    return EXIT_SUCCESS;
}

void get_kmer_totals( hash_table_t *target_kmers,
                      sequence_t **designed_oligos,
                      int num_oligos, int num_mismatches
                    )
{
    hash_table_t *target_copy  = NULL;
    hash_table_t *target_ptr   = target_kmers;
    HT_Entry **items           = NULL;
    hash_table_t *subset_kmers = NULL;
    char       *current_oligo  = NULL;


    unsigned int index = 0;
    
    items = ht_get_items( target_ptr );

    #pragma omp parallel shared( target_ptr, items, designed_oligos ) \
            private( index, target_copy, current_oligo, subset_kmers )
    {
        int oligo_size = designed_oligos[ 0 ]->sequence->size; // assume all oligos are same size
        int num_subsets = num_substrings( oligo_size, WINDOW_SIZE );

        unsigned int inner_index = 0;
        kmer_t *val_copy = NULL;

        HT_Entry **my_items     = NULL;
        HT_Entry **subset_items = NULL;
        HT_Entry *current_item  = NULL;

        int *current_val = NULL;

        target_copy  = malloc( sizeof( hash_table_t ) );
        subset_kmers = malloc( sizeof( hash_table_t ) );

        ht_init( target_copy, LARGE_TABLE_SIZE );

        for( index = 0; index < target_ptr->size; index++ )
            {
                val_copy = malloc( sizeof( kmer_t ) );
                copy_kmer( val_copy, items[ index ]->value );
                ht_add( target_copy, items[ index ]->key,
                        val_copy );
            }

        #pragma omp for
        for( index = 0; index < (unsigned int) num_oligos; index++ )
            {
                ht_init( subset_kmers, num_subsets );
                current_oligo = designed_oligos[ index ]->sequence->data;

                subset_lists_ht( subset_kmers, current_oligo,
                                 oligo_size, WINDOW_SIZE
                               );

                subset_items = ht_get_items( subset_kmers );

                for( inner_index = 0; inner_index < subset_kmers->size; inner_index++ )
                    {


                        get_mismatch_counts( target_copy, items, subset_items[ inner_index ]->key,
                                             target_copy->size, num_mismatches
                                           );
                        free( ( (kmer_t*)(subset_items[ inner_index ]->value) )->seq );
                        free( subset_items[ inner_index ]->value );
                    }
                free( subset_items );
                ht_clear( subset_kmers );
            }

        free( subset_kmers );

        my_items = ht_get_items( target_copy );
        #pragma omp critical
        {
            for( index = 0; index < target_kmers->size; index++ )
                {
                    current_item = my_items[ index ];
                    current_val = (int*) ht_find( target_kmers, current_item->key );
                    *current_val += *(int*) current_item->value;
                }
        }

        for( index = 0; index < target_kmers->size; index++ )
            {
                free( ((kmer_t*)(my_items[ index ]->value))->seq );
                free( my_items[ index ]->value );
            }
        free( my_items );

        ht_clear( target_copy );
        free( target_copy );
    }

    free( subset_kmers );
    free( items );

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

static void subset_lists_local( kmer_t **dest_arr, char *seq,
                                int sequence_len, const int window_size )
{
    int seq_len = sequence_len;
    int num_substr = num_substrings( seq_len, window_size );
    int index = 0;

    unsigned int start      = 0;
    unsigned int end        = 0;
    unsigned const int KMER_SCORE = 0;

    char *substr     = NULL;
    kmer_t *new_kmer = NULL;

    for( index = 0; index < num_substr; index++ )
        {
            substr   = malloc( sizeof( char ) * window_size + 1 );
            new_kmer = malloc( sizeof( kmer_t ) );

            start = index;
            end   = index + window_size;

            substr[ 0 ] = '\0';

            substring_indices( seq, substr, start,
                               end
                             );

            kmer_init( new_kmer, substr, start, end, KMER_SCORE );

            dest_arr[ index ] = new_kmer;
        }
}

static void subset_lists_ht( hash_table_t *dest, char *seq,
                              int sequence_len, const int window_size )
{
    int seq_len = sequence_len;
    int num_substr = num_substrings( seq_len, window_size );
    int index = 0;

    unsigned int start = 0;
    unsigned int end   = 0;
    
    char *substr;

    kmer_t *new_kmer = NULL;

    for( index = 0; index < num_substr; index++ )
        {
            substr = malloc( sizeof( char ) * window_size + 1 );
            substr[ 0 ] = '\0';

            start = index;
            end   = index + window_size;
            substring_indices( seq, substr, start,
                               end
                             );

            if( !ht_find( dest, substr ) )
                {
                    new_kmer = malloc( sizeof( kmer_t ) );

                    kmer_init( new_kmer, substr,
                               start, end, 0
                             );
                    ht_add( dest, new_kmer->seq, new_kmer );
                }
        }
}

hash_table_t *seqs_to_kmer_table( sequence_t **seqs, const int num_seqs )
{
    hash_table_t *table = NULL;
    kmer_t **kmer_arr   = NULL;

    int index       = 0;
    int num_subsets = 0;

    table = malloc( sizeof( hash_table_t ) );
    ht_init( table, LARGE_TABLE_SIZE );

    for( index = 0; index < num_seqs; index++ )
        {
            num_subsets = num_substrings( strlen( seqs[ index ]->sequence->data ), WINDOW_SIZE );
            kmer_arr = malloc( sizeof( kmer_t*) * num_subsets );
            subset_lists_local( kmer_arr, seqs[ index ]->sequence->data,
                                seqs[ index ]->sequence->size, WINDOW_SIZE );

            add_valid_kmers( kmer_arr, num_subsets, table );

            free( kmer_arr );
        }

    return table;
}
void get_mismatch_counts( hash_table_t *table, HT_Entry **items, char *kmer,
                          unsigned int num_items, int num_mismatches
                     )

{
    unsigned int index = 0;
    unsigned int oligo_size = strlen( items[ 0 ]->key );
    kmer_t *value = NULL;
    HT_Entry *current_item = NULL;

    for( index = 0; index < num_items; index++ )
        {
            current_item = items[ index ];

            if( tolerable_match( current_item->key, kmer, oligo_size, num_mismatches ) )
                {
                    value = (kmer_t*) ht_find( table, current_item->key );
                    value->kmer_score++;

                }
            
        }
}

void kmer_init( kmer_t *kmer, char *seq,
                unsigned int start,
                unsigned int end, unsigned int score
              )
{
    kmer->seq = seq;

    kmer->kmer_start = start;
    kmer->kmer_end   = end;
    kmer->kmer_score = score;
}

void write_outputs( char *out_file, hash_table_t *table )
{
    FILE *open_file = fopen( out_file, "w" );
    HT_Entry **ht_items = ht_get_items( table );
    HT_Entry *current_item = NULL;
    kmer_t *current_kmer   = NULL;
    
    unsigned int index = 0;
    const char *HEADER = "Kmer\tScore\tStart\tEnd";

    fprintf( open_file, "%s\n", HEADER );

    for( index = 0; index < table->size; index++ )
        {
            current_item = ht_items[ index ];
            current_kmer = current_item->value;

            fprintf( open_file, "%s\t%u\t%u\t%u\n",
                     current_item->key,
                     current_kmer->kmer_score,
                     current_kmer->kmer_start,
                     current_kmer->kmer_end
                   );
        }




    fclose( open_file );
    free( ht_items );
    
    
}

void clear_table( hash_table_t *table )
{
    HT_Entry **items = NULL;
    unsigned int index = 0;

    items = ht_get_items( table );

    for( index = 0; index < table->size; index++ )
        {
            free( items[ index ]->value );
        }


    ht_clear( table );
    free( table );
    free( items );
}

static inline void add_valid_kmers( kmer_t **kmers, const unsigned int num_subsets, hash_table_t *table )
{
    unsigned int index;
    kmer_t *current_kmer = NULL;

    bool delete_flag = false;

    for( index = 0; index < num_subsets; index++ )
        {
            current_kmer = kmers[ index ];
            if( strchr( current_kmer->seq, 'X' ) == NULL )
                {
                    if( !ht_find( table, current_kmer->seq ) )
                        {
                            ht_add( table,
                                    current_kmer->seq,
                                    current_kmer
                                    );
                        }
                    else
                        {
                            delete_flag = true;
                        }
                }
            else
                {
                    delete_flag = true;
                }

            if( delete_flag )
                {
                    delete_flag = false;

                    free( current_kmer->seq );
                    free( current_kmer );
                }
        }
}


static inline void copy_kmer( kmer_t *dest, kmer_t *src )
{
    dest->seq = malloc( sizeof( char ) * strlen( src->seq ) + 1 );
    strcpy( dest->seq, src->seq );
    dest->kmer_start = src->kmer_start;
    dest->kmer_end   = src->kmer_end;
    dest->kmer_score = src->kmer_score;
}
