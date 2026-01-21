#include "../include/common.h"
#include "../include/barcode_match.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/utils.h"
#include "../include/memory.h"
#include "../include/io.h"
#include "../include/plot_histogram.h"
// EMfit.h removed - EM functionality no longer needed
#include "../include/heatmap.h"

//will  print if DEBUG is set or debug=1
//code for feature sequences stats

//code for heatmap generation
#ifndef NO_HEATMAP
// Heatmap functionality is now in src/heatmap.c
#endif

void destroy_feature_counts(gpointer data) {
    feature_counts *fc = (feature_counts*)data;
    if (fc && fc->counts) {
        g_hash_table_destroy(fc->counts);
    }
    // fc itself is pool-allocated, so we don't free it here.
}

void destroy_feature_umi_counts(gpointer data) {
    feature_umi_counts *umi_counts = (feature_umi_counts*)data;
    if (umi_counts && umi_counts->counts) {
        g_hash_table_destroy(umi_counts->counts);
    }
    // We do NOT free(umi_counts) itself, because it was allocated
    // from our custom memory pool, which is freed all at once later.
}

/* moved to prototypes/common; provide non-static definition here for linkage */
/* moved to barcode_match.c */
void initialize_statistics(statistics *stats) {
    stats->start_time = get_time_in_seconds();
    stats->nMismatches = 0;
    stats->recovered = 0;
    stats->pending = 0;
    stats->valid = 0;
    stats->pending_recovered = 0;
    stats->total_unmatched_features = 0;
    stats->number_of_reads = 0;
    stats->unmatched_list.first_entry = NULL;
    stats->unmatched_list.last_entry = NULL;
}

void free_unmatched_barcodes_features_list(unmatched_barcodes_features_block_list *list) {
    // The blocks are allocated from a memory pool, which is freed all at once.
    // We don't need to free each block individually. Just reset the list pointers.
    list->first_entry = NULL;
    list->last_entry = NULL;
}
double get_time_in_seconds() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + (time.tv_usec / 1000000.0);
}


void lowerCaseDifferences(char *ref, char *query, int length){
    for (int i=0; i<length; i++){
        if (ref[i] != query[i]){
            if (query[i] == '\n'){
               //put x until the length of the query and then add a \0
               for (int j=i; j<length; j++){
                   query[j]='x';
               }
               query[length]='\0';
               return;
            }
            query[i]=tolower(query[i]);
        }
    } 
}   
int compare_feature_sequences(const void *a, const void *b) {
    const feature_sequences *fa = *(const feature_sequences **)a;
    const feature_sequences *fb = *(const feature_sequences **)b;

    // First sort by feature_index
    if (fa->feature_index != fb->feature_index)
        return (fa->feature_index > fb->feature_index) - (fa->feature_index < fb->feature_index);

    // Then sort by hamming_distance
    if (fa->hamming_distance != fb->hamming_distance)
        return fa->hamming_distance - fb->hamming_distance;

    // Lastly sort by counts
    return fb->counts - fa->counts;
}

/* moved to barcode_match.c */


/* moved to barcode_match.c */

void initialize_unit_sizes(){
    // Size of feature_counts (rounded up to alignment of uint16_t)
    size_t feature_counts_size = sizeof(feature_counts);
    size_t feature_counts_alignment = __alignof__(feature_counts);
    dynamic_struct_sizes.feature_counts = (feature_counts_size + feature_counts_alignment - 1) & ~(feature_counts_alignment - 1);


    // Size of feature_umi_counts (rounded up to alignment of uint16_t)
    size_t feature_umi_counts_size = sizeof(feature_umi_counts);
    size_t feature_umi_counts_alignment = __alignof__(feature_umi_counts);
    dynamic_struct_sizes.feature_umi_counts = (feature_umi_counts_size + feature_umi_counts_alignment - 1) & ~(feature_umi_counts_alignment - 1);


    // Size of feature_sequences (rounded up to alignment of char)
    size_t feature_sequences_size = sizeof(feature_sequences) + maximum_feature_length + 1;
    DEBUG_PRINT( "Feature sequences size %ld\n", feature_sequences_size);
    size_t feature_sequences_alignment = __alignof__(char);  // Use __alignof__ instead of alignof
    dynamic_struct_sizes.feature_sequences = (feature_sequences_size + feature_sequences_alignment - 1) & ~(feature_sequences_alignment - 1);
    DEBUG_PRINT( "Adjusted feature sequences size %ld\n", dynamic_struct_sizes.feature_sequences);

    // Size of unmatched_barcodes_features_block (rounded up to alignment of unsigned char)
    size_t unmatched_barcodes_features_block_size = sizeof(unmatched_barcodes_features_block)
            + sizeof(uint32_t) // For feature_index
            + sizeof(uint16_t) // For match_position
            + barcode_code_length 
            + umi_code_length
            + 1 // For number_of_closest_barcodes
            + (max_barcode_mismatches + 1) * barcode_code_length
            + (max_barcode_mismatches + 1);
    size_t unmatched_barcodes_features_block_alignment = __alignof__(uint32_t);  // Use __alignof__ instead of alignof
    dynamic_struct_sizes.unmatched_barcodes_features_block = (unmatched_barcodes_features_block_size + unmatched_barcodes_features_block_alignment - 1) & ~(unmatched_barcodes_features_block_alignment - 1);
}
int is_directory(const char *path) {
    struct stat path_stat;
    
    // Check if stat call is successful
    if (stat(path, &path_stat) != 0) {
        // Error: path not found or other issue
        perror("Error accessing path");
        return -1;
    }
    
    // Check if it is a directory
    if (S_ISDIR(path_stat.st_mode)) {
        return 1;
    }

    // Check if it is a file
    if (S_ISREG(path_stat.st_mode)) {
        return 0;
    }

    // In case it's neither a file nor a directory
    return -1;
}
void read_unmatched_barcodes_features_block(unmatched_barcodes_features_block *entry_block, unmatched_barcodes_features *entry){
    entry->next=entry_block->next;
    unsigned char *storage = entry_block->storage;
    memcpy(&entry->feature_index, storage, sizeof(uint32_t));
    storage += sizeof(uint32_t);
    memcpy(&entry->match_position, storage, sizeof(uint16_t));
    storage += sizeof(uint16_t);
    entry->barcode=storage;
    entry->umi=entry->barcode+barcode_code_length;
    entry->number_of_closest_barcodes=entry->umi[umi_code_length];
    entry->closest_barcodes=entry->umi+umi_code_length+1;
    entry->Qscores=entry->closest_barcodes+(barcode_code_length)*(max_barcode_mismatches+1);
}
int insert_feature_sequence(char *sequence, uint32_t feature_index, unsigned char hamming_distance, uint16_t match_position, data_structures *hashes, memory_pool_collection *pools){
    gpointer *value=g_hash_table_lookup(hashes->unique_features_match, sequence);
    if (value){
        feature_sequences *entry = (feature_sequences*)value;
        entry->counts++;
        return 0;
    }
    else{
        feature_sequences *new_entry = (feature_sequences*) allocate_memory_from_pool(pools->feature_sequences_pool);
        strcpy(new_entry->sequence, sequence);
        new_entry->feature_index=feature_index;
        new_entry->hamming_distance=hamming_distance;
        new_entry->match_position=match_position;
        new_entry->counts=1;
        g_hash_table_insert(hashes->unique_features_match, new_entry->sequence, new_entry);
        return 1;
    }   
}
int print_feature_sequences(feature_arrays *features, int *total_counts, char *directory, data_structures *hashes){
    //remember that 1 is the first feature
    memset(total_counts, 0, features->number_of_features * sizeof(int));
    char feature_sequences_filename[FILENAME_LENGTH];
    sprintf(feature_sequences_filename, "%s/feature_sequences.txt", directory);
    FILE *feature_sequencesfp = fopen(feature_sequences_filename, "w");
    if (feature_sequencesfp == NULL) {
        fprintf(stderr, "Failed to open feature sequences file\n");
        exit(EXIT_FAILURE);
    }   
    GHashTableIter iter;
    gpointer key, value;
    int i=0;
    g_hash_table_iter_init(&iter, hashes->unique_features_match);
    int count = g_hash_table_size(hashes->unique_features_match);
    gpointer *array = g_new(gpointer, count);
    if(!array){
        fprintf(stderr, "Failed to allocate memory for array\n");
        exit(EXIT_FAILURE);
    }   
    while (g_hash_table_iter_next(&iter, &key, &value)) { 
        array[i++] = value;
    }
    qsort(array, count, sizeof(gpointer), compare_feature_sequences);
    fprintf(feature_sequencesfp, "Feature Index Sequence Hamming Distance Counts Match Position Feature Name\n");   
    for (i = 0; i < count; i++) {
        char annotated_sequence[LINE_LENGTH];
        const int mapped_index=((feature_sequences*)array[i])->feature_index-1;
        feature_sequences *entry = (feature_sequences*)array[i];
        total_counts[mapped_index]+=entry->counts;
        strcpy(annotated_sequence, entry->sequence);
        int feature_length=features->feature_lengths[entry->feature_index-1];
        lowerCaseDifferences(features->feature_sequences[entry->feature_index-1],annotated_sequence,feature_length);
        fprintf(feature_sequencesfp, "%7u %s %2d %7d %5u %s\n", entry->feature_index,annotated_sequence, entry->hamming_distance,entry->counts, entry->match_position, features->feature_names[entry->feature_index-1]);
    }
    g_free(array);
    return 0;
}



unmatched_barcodes_features_block* add_unmatched_barcode_store_feature(unsigned char *barcodes, unsigned char* corrected_barcodes, char *umi, unsigned char *qscores, uint32_t feature_index, int number_of_variants, uint16_t match_position, memory_pool_collection *pools, statistics *stats){
    unmatched_barcodes_features_block *new_entry = (unmatched_barcodes_features_block*)allocate_memory_from_pool(pools->unmatched_barcodes_features_block_pool); 

    if (stats->unmatched_list.first_entry== NULL){
        stats->unmatched_list.first_entry=new_entry;
        stats->unmatched_list.last_entry=new_entry;
    }
    else{
        stats->unmatched_list.last_entry->next=new_entry;
        stats->unmatched_list.last_entry=new_entry;
    }
    unsigned char *storage = new_entry->storage; 
    memcpy(storage, &feature_index, sizeof(uint32_t));
    storage += sizeof(uint32_t);
    memcpy(storage, &match_position, sizeof(uint16_t));
    storage += sizeof(uint16_t);
    memcpy(storage, barcodes, barcode_code_length);
    storage+=barcode_code_length;
    string2code(umi, umi_length, storage);
    storage+=umi_code_length;    
    storage[0]=number_of_variants;
    storage++;
    //leave empty space for the closest barcodes    
    unsigned char* qscores_storage = storage + (barcode_code_length)*(max_barcode_mismatches+1);
    memcpy(storage, corrected_barcodes, number_of_variants * barcode_code_length);
    memcpy(qscores_storage, qscores, number_of_variants);
    new_entry->next=NULL;
    return new_entry;
}


unsigned char* read_whiteList(char *whitelist_filename,GHashTable *hash, int reverse_complement_flag){
    FILE *whitelist_file = fopen(whitelist_filename, "r");
    if (whitelist_file == NULL) {
        fprintf(stderr, "Failed to open whitelist file %s\n", whitelist_filename);
        exit(EXIT_FAILURE);
    }
    char line[LINE_LENGTH];
    //just in case the file is corrupted - calculate the sizes by counting lines
    //find the number of lines in the file by counting the number of newlines
    char ch=0,lastChar=0;
    size_t line_count=0;
    //from the first line calculate the barcode length and check if the sequence is valid
    int sequence_count=0;
    while((ch = fgetc(whitelist_file)) != EOF) {
        if(ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N'){
            sequence_count++;
        }
        else if (ch == '\n') {
            if(sequence_count){
                if (barcode_length==0){
                    barcode_length=sequence_count;
                    barcode_code_length=(barcode_length+3)/4;
                }
                else if (sequence_count != barcode_length){
                    fprintf(stderr, "Error: Invalid barcode length %d\n", sequence_count);
                    exit(EXIT_FAILURE);
                }
                sequence_count=0;
                line_count++;
            }
        }
        else{
            fprintf(stderr, "Error: Invalid character %c in whitelist file\n", ch);
            exit(EXIT_FAILURE);
        }
    }
    // If the last character isn't a newline, increment the count
    if (lastChar != '\n') {
        line_count++;
    }
    
    whitelist = malloc(line_count * barcode_code_length);
    if (whitelist == NULL) {
        fprintf(stderr, "Failed to allocate memory for whitelist storage\n");
        exit(EXIT_FAILURE);
    }
    memset(whitelist, 0, line_count * barcode_code_length);
    //reset the file pointer to the beginning of the file
    fseek(whitelist_file, 0, SEEK_SET);
    size_t nBarcodes=0;
    //set line to zero
    memset(line, 0, LINE_LENGTH);
    char barcode_sequence[LINE_LENGTH];
    while ( fgets(line, LINE_LENGTH, whitelist_file) != NULL) {
        if (!check_sequence(line, barcode_length)){
            fprintf(stderr, "Error: Invalid barcode sequence of expected length %d %s\n",barcode_length,line);
            exit(EXIT_FAILURE);
        }
        //up to 16 characters will be stored in uint32_t
        memset(line+barcode_length, 0, LINE_LENGTH-barcode_length);
        if (reverse_complement_flag){
            reverse_complement_sequence(line, barcode_sequence, barcode_length);
        }
        else{
            memcpy(barcode_sequence, line, barcode_length+1);
        }
        int j=0;
        int i=0;
        const size_t offset=nBarcodes*barcode_code_length;
        memset(whitelist+offset, 0, barcode_code_length);
        while(j<barcode_length){
            unsigned char char_value=seq2code[(unsigned char)barcode_sequence[j]]<<6 | seq2code[(unsigned char)barcode_sequence[j+1]]<<4 | seq2code[(unsigned char)barcode_sequence[j+2]]<<2 | seq2code[(unsigned char)barcode_sequence[j+3]];
            whitelist[offset+i]=char_value;
            i++;
            j+=4;
        }
        if(!g_hash_table_insert(hash,  (uint32_t*)(whitelist+offset) , whitelist+offset)){
            fprintf(stderr, "Error: Failed to insert barcode %s into the hash table %ld\n", line, nBarcodes);
            exit(EXIT_FAILURE); 
        }
        nBarcodes++;
    }
    fprintf(stderr, "Read %ld barcodes\n", nBarcodes);
    fclose(whitelist_file);
    return whitelist;
}
/* moved to barcode_match.c */


/* moved to barcode_match.c */
/* moved to barcode_match.c */
void string2all_codes(char *string, unsigned char codes[][LINE_LENGTH/2+1], int *lengths){
    //4 codes are returned for the string to capture all the possible frames
    char offset[3][LINE_LENGTH];
    const int seqlength = strlen(string);

    strcpy(offset[0], string+1);
    strcpy(offset[1], string+2);
    strcpy(offset[2], string+3);
    lengths[0]=string2code(string, seqlength, codes[0]);

    for (int i=0; i<3; i++){
        lengths[i+1]=string2code(offset[i], seqlength-1-i, codes[i+1]);
    }   
}

int check_barcode(char *barcode, unsigned char *code){
    for (int i=0; i<barcode_code_length; i++){
        string2code(barcode, barcode_length, code);
    }
    if (g_hash_table_lookup(whitelist_hash,code) == NULL){
        return 1;
    }
    return 0;
}
int find_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t feature_code_length, size_t feature_length){
    unsigned char a1[feature_code_length];
    unsigned char *b=feature_code;
    memcpy(a1, sequence_code, feature_code_length);
    unsigned char right_overhang=feature_length % 4;

    if (right_overhang){
        a1[feature_code_length-1] = a1[feature_code_length-1] & (0xff << (8-2*right_overhang));        
    }
    int length=feature_code_length;
    unsigned char *a=a1;
    if(length >= 8){
        uint64_t *a64=(uint64_t*) a;
        uint64_t *b64=(uint64_t*) b;
        for (int i=0; i<length/8; i++){
            if (a64[i] != b64[i]){
                return feature_code_length+1;
            }
        }
        if (length % 8 == 0){
            return 0;
        }
        a+=8*(length/8);
        b+=8*(length/8);
        length=length%8;
    }
    if(length >= 4){
        const uint32_t *a32=(uint32_t*) a;
        const uint32_t *b32=(uint32_t*) b;
        for (int i=0; i<length/4; i++){
            if (a32[i] != b32[i]){
                return feature_code_length+1;
            }
        }
        if (length % 4 == 0){
            return 0;
        }
        a+=4*(length/4);
        b+=4*(length/4);
        length=length%4;
    }
    if (length >= 2){
        const uint16_t *a16=(uint16_t*) a;
        const uint16_t *b16=(uint16_t*) b;
        for (int i=0;i<length/2; i++){
            if (a16[i] != b16[i]){
                return feature_code_length+1;
            }
        }
        if (length % 2 == 0){
            return 0;
        }
        a+=2*(length/2);
        b+=2*(length/2);
        length=length%2;
    }
    for (int i=0;i<length; i++){
        if (a[i] != b[i]){
            return length+1;
        }
    }
    return 0;
}
int fuzzy_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t feature_code_length, size_t feature_length, int maxHammingDistance){
    int mismatches=0;
    unsigned char a1[feature_code_length];
    unsigned char *b=feature_code;
    memcpy(a1, sequence_code, feature_code_length);
    unsigned char right_overhang=feature_length % 4;
    if (right_overhang){
        a1[feature_code_length-1] = a1[feature_code_length-1] & (0xff << (8-2*right_overhang));        
    }
    unsigned char *a=a1;
    int length=feature_code_length;
    if(length >= 8){
        const uint64_t *a64=(uint64_t*) a;
        const uint64_t *b64=(uint64_t*) b;
        for (int i=0; i<length/8; i++){
            if (a64[i] != b64[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint64_t diff=a64[i]^b64[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]]+diff2Hamming[diff8[2]]+diff2Hamming[diff8[3]] + diff2Hamming[diff8[4]]+diff2Hamming[diff8[5]]+diff2Hamming[diff8[6]]+diff2Hamming[diff8[7]];
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }
            }
        }
        if (length % 8 == 0){
            return mismatches;
        }
        a+=8*(length/8);
        b+=8*(length/8);
        length=length%8;
    }
    if (length >= 4 ){
        const int last=8*(length/8);
        const uint32_t *a32=(uint32_t*) a;
        const uint32_t *b32=(uint32_t*) b;
        for (int i=last/4; i<length/4; i++){
            if (a32[i] != b32[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint32_t diff=a32[i]^b32[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]]+diff2Hamming[diff8[2]]+diff2Hamming[diff8[3]];
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }
            }
        }
        if (length % 4 == 0){
            return mismatches;
        }
        a+=4*(length/4);
        b+=4*(length/4);
        length=length%4;
        
    }
    if (length >= 2){
        const uint16_t *a16=(uint16_t*) a;
        const uint16_t *b16=(uint16_t*) b;
        const int last=4*(length/4);
        for (int i=last/2;i<length/2; i++){
            if (a16[i] != b16[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint16_t diff=a16[i]^b16[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]];  
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }   
            }
        }
        if (length % 2 == 0){
            return mismatches;
        }
        a+=2*(length/2);
        b+=2*(length/2);
        length=length%2;
    }
    const int last=2*(length/2);
    for (int i=last;i<length; i++){
        if (a[i] != b[i]){
            if(mismatches==maxHammingDistance){
                return feature_code_length+1;
            }
            mismatches+=diff2Hamming[a[i]^b[i]];
            if (mismatches > maxHammingDistance){
                return feature_code_length+1;
            }
        }
    }   
    return mismatches;
}

int find_matches_in_sub_arrays(unsigned char *sequence_code, unsigned char *feature_code, size_t sequence_code_length, size_t feature_code_length, size_t feature_length, int maxHammingDistance, int *best_offset){
    int offset=0;
    int length= (sequence_code_length > feature_code_length) ? sequence_code_length : feature_code_length;
    int minHammingDistance=length+1;
    *best_offset=0;
    while (sequence_code_length - offset >= feature_code_length){
        if (!maxHammingDistance){
            if (!find_matching_codes(sequence_code+offset, feature_code, feature_code_length, feature_length)){
                *best_offset=offset;
                return 0;
            }
        }
        else{
            int hammingDistance=fuzzy_matching_codes(sequence_code+offset, feature_code, feature_code_length, feature_length, maxHammingDistance);
            if (hammingDistance < minHammingDistance){
                minHammingDistance=hammingDistance;
                *best_offset=offset;
            }
            //if we have found a perfect match return immediately
            if  (!minHammingDistance){
                *best_offset=offset;
                return 0;
            }   
        }
        offset++;
    }
    return minHammingDistance;
}

int checkSequenceAndCorrectForN(char *line, char *corrected_lines[], char *buffer,int sequence_length, int maxN){
    int nCount=0;
    int indices[maxN];
    corrected_lines[0]=line;
    for (int i=0; i<sequence_length; i++){
        if (line[i] == 'N'){
            indices[nCount++]=i;
            if (nCount >=  maxN) return 0;
            
        }
        else if (line[i] != 'A' && line[i] != 'C' && line[i] != 'G' && line[i] != 'T'){
            return 0;
        }
    }
    if (!nCount){
        return 1;
    }
    char line_copy[LINE_LENGTH];
    strcpy(line_copy, line);
    return generate_sequences(line_copy, sequence_length, indices, &buffer, 0, nCount, corrected_lines, 0);
}

char* printToBuffer(char *string, int sequence_length, char *buffer){
    for (int i=0; i<sequence_length; i++){
        buffer[i]=string[i];
    }
    buffer[sequence_length]='\0';
    return buffer + sequence_length+1;
}

int generate_sequences(char* string, int string_length, int *indices, char **output, int index, int number_of_indices, char *output_indices[], 
int number_of_outputs) {
    static char bases[] = "ACGT";
    //base case when index is equal to number of indices
    if (index == number_of_indices) {
        printToBuffer(string, string_length, *output);
        output_indices[number_of_outputs++]=*output;
        *output += string_length + 1;
        return number_of_outputs;
    }
    //when the index is less than the number of indices
    for (int i = 0; i < 4; i++) {
        string[indices[index]] = bases[i];  
        number_of_outputs=generate_sequences(string, string_length, indices, output, index+1, number_of_indices, output_indices, number_of_outputs);
    }
    return number_of_outputs;
}
int simple_code_search(feature_arrays *features, char *sequence){
    const int string_length=strlen(sequence);
    unsigned char code[string_length/2+1];
    memset(code, 0, string_length/2+1);
    const int sequence_code_length=string2code(sequence,string_length, code);
    for (int i=0; i<features->number_of_features; i++){
        const int feature_code_length=features->feature_code_lengths[i];
        if (sequence_code_length != feature_code_length){
            continue;
        }
        if (find_matching_codes(code, features->feature_codes[i], sequence_code_length, features->feature_lengths[i])== 0){
            return i+1;
        }
    }
    return 0;
}
int simple_hash_search(feature_arrays *features, char *sequence){
    const int string_length=strlen(sequence);
    const int feature_sequence_length=features->feature_lengths[0];
    const int sequence_code_length=features->feature_code_lengths[0];
    if (string_length < feature_sequence_length) return 0;
    char local_line[LINE_LENGTH];
    memcpy(local_line, sequence, feature_sequence_length);
    local_line[feature_sequence_length]='\0';
    unsigned char code[sequence_code_length];
    string2code(local_line, feature_sequence_length, code);
    
    // Create a GBytes key for lookup.
    GBytes *key = g_bytes_new(code, sequence_code_length);
    gpointer value = g_hash_table_lookup(feature_code_hash, key);
    g_bytes_unref(key); // Free the GBytes wrapper

    if (value){
        // The value is the feature index, convert it back from a pointer.
        return GPOINTER_TO_INT(value);
    }
    //check against the mismatched features
    for (int i=0; i<features->number_of_mismatched_features; i++){
        char *query=sequence;
        char *feature=features->feature_sequences[features->mismatched_feature_indices[i]];
        while (*query == *feature){
            query++;
            feature++;
        }
        if (!*feature){
            return features->mismatched_feature_indices[i]+1;
        }
    }
    return 0;
}

int simple_search(feature_arrays *features, char *line){
    if (features->number_of_features < 150){
        for (int i=0; i<features->number_of_features; i++){
            char *query=line;
            char *feature=features->feature_sequences[i];
            //DEBUG_PRINT( "Comparing %s %s\n", query, feature);
            while (*query == *feature){
                query++;
                feature++;
            }
            if (!*feature){
                return i+1;
            }
        }
        return 0;
    }
    return simple_hash_search(features, line);    
}

int simple_hamming_search(feature_arrays *features, char *line, int maxHammingDistance, int *hamming_distance){  
    int ambiguous=0;
    int bestFeature=0;
    int bestHammingDistance=maxHammingDistance;
    for (int i=0; i<features->number_of_features; i++){
        char *query=line;
        char *feature=features->feature_sequences[i];
        int hammingDistance=0;
        
        //DEBUG_PRINT( "Comparing %s %s\n", query, feature);
        while (*query && *feature){
            if (*query != *feature){
                hammingDistance++;
                if (hammingDistance > maxHammingDistance){
                    break;
                }
            }
            query++;
            feature++;
        }
        if (!*feature){
            if (hammingDistance == bestHammingDistance){
                if(bestFeature){
                    ambiguous=1;
                }
                else{
                    bestFeature=i+1;
                    bestHammingDistance=hammingDistance;
                    ambiguous=0;
                }
            }
            else if (hammingDistance < bestHammingDistance){
               bestFeature=i+1;
               bestHammingDistance=hammingDistance;
               ambiguous=0;
            }
        }

    }
    *hamming_distance=bestHammingDistance;
    if (bestHammingDistance <= maxHammingDistance && !ambiguous){
        return bestFeature;
    }
    return 0;
}
int find_feature_match_single(feature_arrays *features, char *lineR2, int maxHammingDistance,int *bestScore, char **matching_sequence, uint16_t *match_position){
    // convert lineR2 to 4 codes
    // do a quick check to see if there is a perfect match the constant feature
    //int best_feature=simple_search_code(features, lineR2, constant_offset);
    unsigned char codes[4][LINE_LENGTH/2+1];
    int code_lengths[4];
    string2all_codes(lineR2, codes, code_lengths);
    int best_feature=0;
    int bestHammingDistance=maxHammingDistance; 
    int best_sequence_offset=0;
    int ambiguous=0;
    for (int i=0; i<4; i++){
        for (int j=0; j<features->number_of_features; j++){
            int code_offset=0;
            int hammingDistance=find_matches_in_sub_arrays(codes[i], features->feature_codes[j], code_lengths[i], features->feature_code_lengths[j],features->feature_lengths[j], maxHammingDistance, &code_offset);
            if (!hammingDistance){
                *bestScore=0;
                best_feature=j+1;
                *matching_sequence=lineR2+i+code_offset*4;
                *match_position = i+code_offset*4;
                return best_feature;
            }
            if (hammingDistance < bestHammingDistance){
                best_feature=j+1;
                bestHammingDistance=hammingDistance;
                best_sequence_offset=i+4*code_offset;
                ambiguous=0;
            }
            else if (hammingDistance == bestHammingDistance){
                ambiguous=1;
            }
        }
    }
    if (bestHammingDistance <= maxHammingDistance){
        *bestScore=bestHammingDistance;
        *matching_sequence=lineR2+best_sequence_offset;
        if (!ambiguous){
            return best_feature;
        }
        return 0;
    }
    return 0;
}
int find_feature_match_parallel(feature_arrays *features, char *lineR2, int maxHammingDistance, int nThreads, int *bestScore, char **matching_sequence, uint16_t *match_position){
    // convert lineR2 to 4 codes
    // do a quick check to see if there is a perfect match the constant feature
    //int best_feature=simple_search_code(features, lineR2, constant_offset);
    nThreads=(nThreads > 4) ? 4 : nThreads;
    nThreads=(nThreads < 1) ? 1 : nThreads;
    if (nThreads==1){
        return find_feature_match_single(features, lineR2, maxHammingDistance, bestScore, matching_sequence, match_position);
    }
    unsigned char codes[4][LINE_LENGTH/2+1];
    int code_lengths[4];
    string2all_codes(lineR2, codes, code_lengths);

    int best_feature=0;
    int bestFeatureDistance=maxHammingDistance;   
    int bestHammingDistances[4]={maxHammingDistance,maxHammingDistance,maxHammingDistance,maxHammingDistance};
    int best_code_offsets[4]={0,0,0,0};
    int ambiguous[4]={0,0,0,0};
    int best_match[4]={0,0,0,0};
    int best_matching_indices[4]={0,0,0,0};

    int exact_match_found=0;
    #pragma omp parallel for num_threads(nThreads)
    for (int i=0; i<4; i++){
        //get the thread number
        const int thread_num=omp_get_thread_num();
        for (int j=0; j<features->number_of_features && !exact_match_found; j++){
            int code_offset=0;
            int hammingDistance=find_matches_in_sub_arrays(codes[i], features->feature_codes[j], code_lengths[i], features->feature_code_lengths[j],features->feature_lengths[j], maxHammingDistance, &code_offset);
            if (!hammingDistance){
                *bestScore=0;
                exact_match_found=1;
                best_feature=j+1;
                *matching_sequence=lineR2+i+code_offset*4;
                *match_position = i+code_offset*4;
                best_matching_indices[thread_num]=i;
                break;
            }
            if (hammingDistance < bestHammingDistances[thread_num]){
                best_match[thread_num]=j+1;
                bestHammingDistances[thread_num]=hammingDistance;
                best_code_offsets[thread_num]=code_offset;
                ambiguous[thread_num]=0;
                best_matching_indices[thread_num]=i;
            }
            else if (hammingDistance == bestHammingDistances[thread_num]){
                if (best_match[thread_num]){
                    ambiguous[thread_num]=1;
                }
                else{
                    ambiguous[thread_num]=0;
                    best_match[thread_num]=j+1;
                    best_code_offsets[thread_num]=code_offset;
                    bestHammingDistances[thread_num]=hammingDistance;
                    best_matching_indices[thread_num]=i;
                }
            }
        }
    }
    if (exact_match_found){//matching_sequence is already set
        return best_feature;
    }
    //find the best hamming distance - if there is a perfect match return immediately
    char multiAmbiguous=0;
    int best_code_offset=0;
    int best_i=0;
    for (int i=0; i<nThreads; i++){
        if (bestHammingDistances[i] < bestFeatureDistance){
            bestFeatureDistance=bestHammingDistances[i];
            best_feature=best_match[i];
            best_code_offset=best_code_offsets[i];
            multiAmbiguous=ambiguous[i]; //reset ambiguous flag
            best_i=best_matching_indices[i];
        }
        else if (bestHammingDistances[i] == bestFeatureDistance){
            if (best_feature ){
                multiAmbiguous=1;
            }
            else{
                best_feature=best_match[i];
                bestFeatureDistance=bestHammingDistances[i];
                best_code_offset=best_code_offsets[i];
                best_i=best_matching_indices[i];
                multiAmbiguous=ambiguous[i];
            }
        }
    }
    if (bestFeatureDistance <= maxHammingDistance){
        *bestScore=bestFeatureDistance;
        *matching_sequence=lineR2+best_i+best_code_offset*4;
        *match_position=best_i+best_code_offset*4;
        if (!multiAmbiguous){
            return best_feature;
        }
        return 0;
    }
    return 0;
}
void update_feature_counts(char *barcodeString, char *umi, uint32_t feature_index, data_structures *hashes,  memory_pool_collection *pools){ 
    unsigned char code[barcode_code_length];
    string2code(barcodeString, barcode_length, code);
    update_feature_counts_from_code(code, umi, feature_index, hashes, pools);
}
// In src/assignBarcodes.c

void update_feature_counts_from_code(unsigned char *code, char *umi, uint32_t feature_index, data_structures *hashes, memory_pool_collection *pools){
    update_umi_counts(code, umi, feature_index, hashes, pools);
    feature_counts *s = g_hash_table_lookup(hashes->filtered_hash, code);

    if (s == NULL) {
        // Logic for a NEW barcode
        s = (feature_counts*) allocate_memory_from_pool(pools->feature_counts_pool);
        if (s == NULL) return; // Error check

        memcpy(s->sequence_code, code, barcode_code_length);

        // NEW: Initialize the inner hash table for this barcode's counts.
        // We don't provide destroy functions because the keys/values are not allocated pointers.
        s->counts = g_hash_table_new(g_direct_hash, g_direct_equal);
        if (s->counts == NULL) {
            fprintf(stderr, "Fatal: GHashTable creation for counts failed.\n");
            exit(EXIT_FAILURE);
        }

        g_hash_table_insert(hashes->filtered_hash, s->sequence_code, s);

        // For a new entry, the count for this feature is 1.
        g_hash_table_insert(s->counts, GUINT_TO_POINTER(feature_index), GUINT_TO_POINTER(1));
        // The total count (at index 0) is also 1.
        g_hash_table_insert(s->counts, GUINT_TO_POINTER(0), GUINT_TO_POINTER(1));

    } else {
        // Logic for an EXISTING barcode
        // 1. Increment the specific feature's count
        uintptr_t current_count = GPOINTER_TO_UINT(g_hash_table_lookup(s->counts, GUINT_TO_POINTER(feature_index)));
        current_count++;
        g_hash_table_replace(s->counts, GUINT_TO_POINTER(feature_index), GUINT_TO_POINTER(current_count));

        // 2. Increment the total count stored at key '0'
        uintptr_t total = GPOINTER_TO_UINT(g_hash_table_lookup(s->counts, GUINT_TO_POINTER(0)));
        total++;
        g_hash_table_replace(s->counts, GUINT_TO_POINTER(0), GUINT_TO_POINTER(total));
    }
}
void update_umi_counts(unsigned char *code, char *umi,  uint32_t feature_index,data_structures *hashes, memory_pool_collection *pools){
    unsigned char code8[8];
    memset(code8,0,8);
    memcpy(code8,code,barcode_code_length);
    string2code(umi, umi_length, code8+barcode_code_length);

    feature_umi_counts *s = g_hash_table_lookup(hashes->sequence_umi_hash, code8);
    
    if (s == NULL) {
        // Logic for a NEW barcode-UMI combination
        s = (feature_umi_counts*) allocate_memory_from_pool(pools->feature_umi_counts_pool);
        if (s == NULL) return; // Error check

        memcpy(s->sequence_umi_code, code8, 8);
        
        // Create the new inner hash table for this UMI's feature counts
        s->counts = g_hash_table_new(g_direct_hash, g_direct_equal);
        
        // Set initial count for this feature to 1
        g_hash_table_insert(s->counts, GUINT_TO_POINTER(feature_index), GUINT_TO_POINTER(1));
        // Mark as unvisited for the connected components algorithm by setting key '0' to 0.
        g_hash_table_insert(s->counts, GUINT_TO_POINTER(0), GUINT_TO_POINTER(0));

        g_hash_table_insert(hashes->sequence_umi_hash, s->sequence_umi_code, s);
    } else {
        // Logic for an EXISTING barcode-UMI combination
        uintptr_t current_count = GPOINTER_TO_UINT(g_hash_table_lookup(s->counts, GUINT_TO_POINTER(feature_index)));
        current_count++;
        g_hash_table_replace(s->counts, GUINT_TO_POINTER(feature_index), GUINT_TO_POINTER(current_count));
    }
}
char check_neighbor(uint64_t code64,uint32_t *counts, data_structures *hashes){
    feature_umi_counts *result = g_hash_table_lookup(hashes->sequence_umi_hash, &code64);

    if (result && result->counts) {
        // Check if this node has been visited using key '0'
        uintptr_t visited_flag = GPOINTER_TO_UINT(g_hash_table_lookup(result->counts, GUINT_TO_POINTER(0)));
        
        if (visited_flag == 0) { // If not visited
            // --- NEW: Iterate through the hash table to sum counts ---
            GHashTableIter iter;
            gpointer key, value;
            g_hash_table_iter_init(&iter, result->counts);
            while (g_hash_table_iter_next(&iter, &key, &value)) {
                uint32_t feature_index = GPOINTER_TO_UINT(key);
                if (feature_index > 0) { // Don't add the visited flag to the temp counts array
                    counts[feature_index] += GPOINTER_TO_UINT(value);
                }
            }

            // Mark as visited by replacing the value at key '0' with 1
            g_hash_table_replace(result->counts, GUINT_TO_POINTER(0), GUINT_TO_POINTER(1));
            return 1; // Indicate success
        }
    }
    return 0; // Not found or already visited
}
int find_neighbors(uint64_t key64, uint64_t *neighbors, uint32_t *counts, data_structures *hashes){
    int neighbor_count=0;
    uint64_t code64=key64;
    unsigned char *code8= (unsigned char*) &code64;
    unsigned char *code=code8+barcode_code_length;
    for (int i=0; i<umi_code_length; i++){       
        const unsigned char mask=code[i];
        for (unsigned char j=1; j<4; j++){
            code[i]=mask ^ j;
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask ^ (j << 2);
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask ^ (j << 4);
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask;
        }
    }
    return neighbor_count;
}
void find_deduped_counts(data_structures *hashes, GHashTable* barcode_to_deduped_counts, uint16_t stringency, uint16_t min_counts){
    GHashTableIter iter;
    gpointer lookup_key, result;
    g_hash_table_iter_init(&iter, hashes->sequence_umi_hash);
    uint32_t clique_counts[number_of_features+1];
    
    // Create tracking set for barcodes that might need cleanup (empty after deduping)
    GHashTable *empty_candidates = g_hash_table_new(g_direct_hash, g_direct_equal);

    while (g_hash_table_iter_next(&iter, &lookup_key, &result)) {
        feature_umi_counts *umi_counts = (feature_umi_counts*) result;
        if (GPOINTER_TO_UINT(g_hash_table_lookup(umi_counts->counts, GUINT_TO_POINTER(0)))) {
            continue; // Already processed if the flag is 1
        }
        memset(clique_counts, 0, sizeof(clique_counts));
        
        find_connected_component(lookup_key, clique_counts, hashes);

        feature_counts *s = g_hash_table_lookup(hashes->filtered_hash, umi_counts->sequence_umi_code);
        if (s) {
            // Get or create the inner hash table for this barcode's deduped counts
            GHashTable* temp_deduped_hash = g_hash_table_lookup(barcode_to_deduped_counts, s->sequence_code);
            if (temp_deduped_hash == NULL) {
                temp_deduped_hash = g_hash_table_new(g_direct_hash, g_direct_equal);
                g_hash_table_insert(barcode_to_deduped_counts, s->sequence_code, temp_deduped_hash);
            }
            
            add_deduped_count(temp_deduped_hash, clique_counts, stringency, min_counts);
            
            // If still empty after deduping, remember this barcode for cleanup
            if (g_hash_table_size(temp_deduped_hash) == 0) {
                g_hash_table_replace(empty_candidates,
                                     s->sequence_code,   /* key */
                                     NULL);              /* value */
            }
        }
    }
    
    // --- FINAL CLEAN-UP SWEEP ---
    // Remove barcodes that are truly empty after all deduping is complete
    GHashTableIter del_iter;
    gpointer bcode_key, unused_val;
    g_hash_table_iter_init(&del_iter, empty_candidates);

    while (g_hash_table_iter_next(&del_iter, &bcode_key, &unused_val)) {
        GHashTable *dedup_hash =
            g_hash_table_lookup(barcode_to_deduped_counts, bcode_key);

        if (dedup_hash && g_hash_table_size(dedup_hash) == 0) {
            /* Remove from main table and free inner hash */
            g_hash_table_steal(barcode_to_deduped_counts, bcode_key);
            g_hash_table_destroy(dedup_hash);
        }
    }
    g_hash_table_destroy(empty_candidates);
}

void find_connected_component(gpointer start_key, uint32_t *counts, data_structures *hashes){
    uint64_t neighbors[umi_length*3];
    clear_queue(hashes->neighbors_queue);
    enqueue(hashes->neighbors_queue,*(uint64_t*)start_key);
    while (!is_empty(hashes->neighbors_queue)) {
        uint64_t code = dequeue(hashes->neighbors_queue);
        check_neighbor(code,counts,hashes);
        int neighbor_count=find_neighbors(code,neighbors, counts, hashes);
        //if (neighbor_count) DEBUG_PRINT( "Found %d neighbors\n",neighbor_count);
        for (int i = 0; i < neighbor_count; i++) {
            enqueue(hashes->neighbors_queue,neighbors[i]);
        }
    }
}
/**
 * @brief Calculates the winning feature based on UMI clique counts and stringency,
 * and updates a temporary hash table with the deduplicated count.
 *
 * @param temp_deduped_hash A temporary hash table mapping a feature_index to its deduped count.
 * @param clique_counts An array of raw UMI counts for the connected component.
 * @param stringency The user-defined stringency for deduplication.
 * @param min_counts The minimum count threshold.
 */
 void add_deduped_count(GHashTable* temp_deduped_hash, uint32_t *clique_counts, uint16_t stringency, uint16_t min_counts) {
    uint32_t winning_feature_index = 0;

    // --- STRINGENCY LOGIC (Restored from original code) ---
    if (!stringency) {
        // RNA-seq strategy: any feature with enough counts gets a single deduped count.
        // This is a special case as it can increment multiple features.
        for (int i = 1; i < number_of_features + 1; i++) {
            if (clique_counts[i] > min_counts) {
                uintptr_t count = GPOINTER_TO_UINT(g_hash_table_lookup(temp_deduped_hash, GUINT_TO_POINTER(i)));
                count++;
                g_hash_table_replace(temp_deduped_hash, GUINT_TO_POINTER(i), GUINT_TO_POINTER(count));
            }
        }
        return; // Return early as the logic is different
    }

    if (stringency >= 1000) {
        // Highest stringency: only one feature can have counts.
        uint32_t feature_index = 0;
        uint32_t total_counts = 0;
        for (int i = 1; i < number_of_features + 1; i++) {
            if (clique_counts[i]) {
                if (feature_index) { // If a feature has already been found, there's more than one.
                    winning_feature_index = 0; // Invalidate winner
                    break;
                }
                feature_index = i;
                total_counts = clique_counts[i];
            }
        }
        if (feature_index && total_counts > min_counts) {
            winning_feature_index = feature_index;
        }
    } else {
        // Mid-stringency: find the feature with the highest, unique count.
        uint32_t feature_index = 0;
        uint32_t max_counts = 0;
        uint32_t total_counts = 0;
        unsigned char unique = 0;
        for (int i = 1; i < number_of_features + 1; i++) {
            if (clique_counts[i]) {
                total_counts += clique_counts[i];
                if (clique_counts[i] > max_counts) {
                    max_counts = clique_counts[i];
                    feature_index = i;
                    unique = 1;
                } else if (clique_counts[i] == max_counts) {
                    unique = 0;
                }
            }
        }
        if (unique && (double)max_counts > total_counts * stringency / 1000.0 && total_counts > min_counts) {
            winning_feature_index = feature_index;
        }
    }

    // --- If a single winner was found, increment its count in the temporary hash table ---
    if (winning_feature_index > 0) {
        uintptr_t current_deduped_count = GPOINTER_TO_UINT(
            g_hash_table_lookup(temp_deduped_hash, GUINT_TO_POINTER(winning_feature_index))
        );
        current_deduped_count++;
        g_hash_table_replace(temp_deduped_hash, 
                             GUINT_TO_POINTER(winning_feature_index), 
                             GUINT_TO_POINTER(current_deduped_count));
    }
}

/* moved to barcode_match.c */

/* Increment (and lazily allocate/resize) the histogram for one feature. */
static void update_feature_hist(GArray **feature_hist,
                                int feature_idx, int count)
{
    GArray *h = feature_hist[feature_idx];
    if (!h)
        h = feature_hist[feature_idx] =
                 g_array_new(FALSE, TRUE, sizeof(uint32_t));

    if (h->len <= (guint)count)
        g_array_set_size(h, count + 1);

    g_array_index(h, uint32_t, count)++;
}

/* Build the cumulative histogram in slot 0 once all features are done. */
static void build_cumulative_feature_hist(feature_arrays *features,
                                          GArray **feature_hist)
{
    feature_hist[0] = g_array_new(FALSE, TRUE, sizeof(uint32_t));
    if (!feature_hist[0])
        return;

    /* find maximum length */
    size_t max_len = 0;
    for (int i = 1; i <= features->number_of_features; ++i)
        if (feature_hist[i] && feature_hist[i]->len > max_len)
            max_len = feature_hist[i]->len;

    g_array_set_size(feature_hist[0], max_len);
    for (guint j = 0; j < max_len; ++j)
        g_array_index(feature_hist[0], uint32_t, j) = 0;

    /* accumulate */
    for (int i = 1; i <= features->number_of_features; ++i) {
        GArray *h = feature_hist[i];
        if (!h) continue;
        for (guint j = 0; j < h->len; ++j)
            g_array_index(feature_hist[0], uint32_t, j) +=
                g_array_index(h, uint32_t, j);
    }
}

void printFeatureCounts(feature_arrays *features, int *deduped_counts, int *barcoded_counts, GArray **feature_hist, char *directory, data_structures *hashes, statistics *stats, GHashTable *barcode_to_deduped_hash, GHashTable *filtered_barcodes_hash){
    int total_deduped_counts = 0;
    int total_raw_counts = 0;
    char barcodes_file[FILENAME_LENGTH];
    char stats_file[FILENAME_LENGTH];

    memset(deduped_counts, 0, (features->number_of_features + 1) * sizeof(int));
    memset(barcoded_counts, 0, (features->number_of_features + 1) * sizeof(int));


    char output_directory[FILENAME_LENGTH];
    strcpy(output_directory, directory);
    if (filtered_barcodes_hash) {
        sprintf(output_directory, "%s/filtered", directory);
    }

    mkdir_p(output_directory);
    sprintf(barcodes_file, "%s/barcodes.txt", output_directory);
    sprintf(stats_file, "%s/stats.txt", output_directory);
    FILE *barcodesfp = fopen(barcodes_file, "w");
    if (barcodesfp == NULL) {
        fprintf(stderr, "Error opening barcodes file %s\n", barcodes_file);
        exit(1);
    }

    char matrix_file[LINE_LENGTH];
    sprintf(matrix_file, "%s/matrix.mtx", output_directory);
    FILE *matrixfp = fopen(matrix_file, "w");
    if (matrixfp == NULL) {
        fprintf(stderr, "Error opening matrix file %s\n", matrix_file);
        exit(1);
    }
    
    // --- Step 2: Calculate stats and write the Matrix Market header ---
    size_t number_of_features_seen = 0;
    size_t number_of_barcode_entries = 0;
    GHashTableIter iter;
    gpointer key, value;

    GArray *feature_richness_hist = g_array_new(FALSE, TRUE, sizeof(uint32_t));

    g_hash_table_iter_init(&iter, barcode_to_deduped_hash);
    while (g_hash_table_iter_next(&iter, &key, &value)) {
        if (filtered_barcodes_hash) {
            char barcode[barcode_length + 1];
            code2string((unsigned char *)key, barcode, barcode_code_length);
            if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
            if (g_hash_table_lookup(filtered_barcodes_hash, barcode) == NULL) {
                continue;
            }
        }
        number_of_barcode_entries++;
        GHashTable *deduped_hash = (GHashTable *)value;
        guint n_features_in_barcode = g_hash_table_size(deduped_hash);
        number_of_features_seen += n_features_in_barcode;

        if (feature_richness_hist && n_features_in_barcode > 0) {
            if (feature_richness_hist->len <= n_features_in_barcode) {
                g_array_set_size(feature_richness_hist, n_features_in_barcode + 1);
            }
            g_array_index(feature_richness_hist, uint32_t, n_features_in_barcode)++;
        }
    }
    
    fprintf(matrixfp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(matrixfp, "%%metadata_json: {\"software_version\": \"assignBarcodes-0.1\", \"format_version\": 1}\n");
    fprintf(matrixfp, "%d %ld %ld\n", features->number_of_features, number_of_barcode_entries, number_of_features_seen);
    
    
    // --- Step 2.5: Calculate co-expression histograms for richness heat-map ---
    // A more direct approach to building the histograms.

    // Step 1: Find max co-occurrence and build an inverted index (feature -> barcodes).
    GHashTable *feature_to_barcodes_hash =
        g_hash_table_new(g_direct_hash, g_direct_equal);
    int global_max_cooccur = 0;

    gpointer       barcode_key, deduped_hash_value;
    int skipped_barcodes = 0;
    int processed_barcodes = 0;
    g_hash_table_iter_init(&iter, barcode_to_deduped_hash);
    while (g_hash_table_iter_next(&iter, &barcode_key, &deduped_hash_value)) {
        char barcode[barcode_length + 1];
        code2string((unsigned char *)barcode_key, barcode, barcode_code_length);
        if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
        if (filtered_barcodes_hash &&
            g_hash_table_lookup(filtered_barcodes_hash, barcode) == NULL){
            skipped_barcodes++;
            continue;
        }
        processed_barcodes++;
        GHashTable *features_in_barcode = (GHashTable *)deduped_hash_value;
        guint       n_feats_in_bc       = g_hash_table_size(features_in_barcode);
        if (n_feats_in_bc == 0) continue;
        if (n_feats_in_bc  > global_max_cooccur) {
            global_max_cooccur = n_feats_in_bc;
        }

        GList *feat_idxs = g_hash_table_get_keys(features_in_barcode);
        for (GList *l = feat_idxs; l; l = l->next) {
            int f_idx = GPOINTER_TO_INT(l->data);
            if (f_idx <= 0) continue;

            GList *barcodes_for_feature =
                g_hash_table_lookup(feature_to_barcodes_hash, l->data);
            barcodes_for_feature = g_list_prepend(barcodes_for_feature, barcode_key);
            g_hash_table_insert(feature_to_barcodes_hash, l->data,
                                barcodes_for_feature);
        }
        g_list_free(feat_idxs);
    }
    fprintf(stderr, "Skipped barcodes: %d\n", skipped_barcodes);
    fprintf(stderr, "Processed barcodes: %d\n", processed_barcodes);
    // Step 2: Allocate and zero the dense histogram matrix.
    int **coexpression_histograms =
        malloc((features->number_of_features + 1) * sizeof(int *));
    for (int i = 0; i <= features->number_of_features; ++i) {
        coexpression_histograms[i] =
                calloc(global_max_cooccur + 1, sizeof(int));
    }

    // Step 3: Populate the histogram matrix using the inverted index.
    GList *active_features = g_hash_table_get_keys(feature_to_barcodes_hash);
    for (GList *l = active_features; l; l = l->next) {
        int   f_idx    = GPOINTER_TO_INT(l->data);
        GList *barcodes = g_hash_table_lookup(feature_to_barcodes_hash, l->data);

        for (GList *b_node = barcodes; b_node; b_node = b_node->next) {
            char *      barcode_str = (char *)b_node->data;
            GHashTable *deduped_hash =
                g_hash_table_lookup(barcode_to_deduped_hash, barcode_str);

            int num_feats          = g_hash_table_size(deduped_hash);

            if (num_feats > 0) {
                coexpression_histograms[f_idx][num_feats]++;
            }
        }
    }
    g_list_free(active_features);

    // Step 4: Set column 0 to the max co-occurrence for each feature.
    for (int i = 1; i <= features->number_of_features; ++i) {
        int row_max = 0;
        for (int j = global_max_cooccur; j >= 1; --j) {
            if (coexpression_histograms[i][j] > 0) {
                row_max = j;
                break;
            }
        }
        coexpression_histograms[i][0] = row_max;
    }

    // --- Clean up ---
    // The barcode strings are owned by barcode_to_deduped_hash, so we only
    // need to free the lists themselves.
    GList *f2b_values = g_hash_table_get_values(feature_to_barcodes_hash);
    for (GList *l = f2b_values; l; l = l->next) {
        g_list_free(l->data);
    }
    g_list_free(f2b_values);
    g_hash_table_destroy(feature_to_barcodes_hash);

    generate_heatmap(output_directory, features, coexpression_histograms);

    for (int i = 0; i <= features->number_of_features; ++i)
        free(coexpression_histograms[i]);
    free(coexpression_histograms);
    
    // --- Step 3: Iterate through all barcodes to write files and calculate final stats ---
    int total_barcodes=0;
    int total_excluded_barcodes=0;
    g_hash_table_iter_init(&iter, hashes->filtered_hash);
    while (g_hash_table_iter_next(&iter, &key, &value)) {
        feature_counts *entry = (feature_counts*)value;
        total_raw_counts += GPOINTER_TO_UINT(g_hash_table_lookup(entry->counts, GUINT_TO_POINTER(0)));
        // --- Populate total_barcoded_counts (raw counts) ---
        GHashTableIter raw_iter;
        gpointer raw_key, raw_value;
        g_hash_table_iter_init(&raw_iter, entry->counts);
        while (g_hash_table_iter_next(&raw_iter, &raw_key, &raw_value)) {
            int feature_index = GPOINTER_TO_INT(raw_key);
            if (feature_index > 0) { // Skip the total stored at index 0
                barcoded_counts[feature_index - 1] += GPOINTER_TO_INT(raw_value);
            }
        }

        // Check if this barcode has any deduplicated counts
        GHashTable* deduped_hash = g_hash_table_lookup(barcode_to_deduped_hash, entry->sequence_code);
        if (deduped_hash && g_hash_table_size(deduped_hash) > 0) {
            
            char barcode[barcode_length + 1];
            code2string(entry->sequence_code, barcode, barcode_code_length);
            if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
            if (filtered_barcodes_hash && g_hash_table_lookup(filtered_barcodes_hash, barcode) == NULL) {
                total_excluded_barcodes++;
                continue;
            }
            
            fprintf(barcodesfp, "%s\n", barcode);
            // Write its matrix entries
            GHashTableIter dedup_iter;
            gpointer dedup_key, dedup_value;
            g_hash_table_iter_init(&dedup_iter, deduped_hash);
            while (g_hash_table_iter_next(&dedup_iter, &dedup_key, &dedup_value)) {
                int deduped_count = GPOINTER_TO_INT(dedup_value);

                int f_idx = GPOINTER_TO_INT(dedup_key);          /* 1-based          */
                deduped_counts[f_idx - 1] += deduped_count;      /* This was missing */
                int c     = deduped_count;                       /* 0,1,2,…          */

                update_feature_hist(feature_hist, f_idx, c);
                fprintf(matrixfp, "%d %d %d\n", GPOINTER_TO_INT(dedup_key), total_barcodes + 1, deduped_count);
                total_deduped_counts += deduped_count;
                
            }
            total_barcodes++;
        }

    }
    
    // --- Step 3.5: Create cumulative histogram at index 0 ---
    build_cumulative_feature_hist(features, feature_hist);
    
    
    
    // --- Step 4: Final Cleanup ---
    // The barcode_to_deduped_hash is managed by the caller (finalize_processing)
    
    fclose(barcodesfp);
    fclose(matrixfp);
    fprintf(stderr,"closing matrix file\n");
    fprintf(stderr,"writing stats file\n");
    FILE *statsfp = fopen(stats_file, "w");
    fprintf (stderr, "Total feature counts %d\n", total_raw_counts);
    fprintf (stderr, "Total deduped feature counts %d\n", total_deduped_counts);
    fprintf (stderr, "Total barcodes %d\n", total_barcodes);
    fprintf (stderr, "Total excluded barcodes %d\n", total_excluded_barcodes);
    fprintf (stderr, "Total unique barcode UMIs %d\n", g_hash_table_size(hashes->sequence_umi_hash));
    fprintf (stderr, "Total whitelisted barcodes %d\n", g_hash_table_size(hashes->filtered_hash));
    fprintf (stderr,"Total feature counts %d total_unmatched_reads %ld\n", total_raw_counts, stats->total_unmatched_features);
    fprintf (stderr, "Percentage reads assigned to barcode %.4f\n", 100.0*(total_raw_counts/(double) (total_raw_counts+stats->total_unmatched_features)));
    fprintf (statsfp, "Total feature counts %d\n", total_raw_counts);
    fprintf (statsfp, "Total deduped feature counts %d\n", total_deduped_counts);
    fprintf (statsfp, "Total barcodes %d\n", total_barcodes);
    fprintf (statsfp, "Total excluded barcodes %d\n", total_excluded_barcodes);
    fprintf (statsfp, "Total unique barcode UMIs %d\n", g_hash_table_size(hashes->sequence_umi_hash));
    fprintf (statsfp, "Total whitelisted barcodes %d\n", g_hash_table_size(hashes->filtered_hash));
    fprintf (statsfp,"Total_unmatched_reads %ld\n", stats->total_unmatched_features);
    fprintf (statsfp, "Percentage reads assigned to barcode %.4f\n", 100.0*(total_raw_counts/(double) (total_raw_counts+stats->total_unmatched_features)));
    fclose(statsfp);

    if (feature_richness_hist) {
        plot_simple_histogram(output_directory, "feature_richness_histogram.html", "Feature Richness per Barcode", "Number of Distinct Features", "Frequency (Barcodes)", feature_richness_hist);
        g_array_unref(feature_richness_hist);
    }
    if (feature_hist[0]) {
        plot_simple_histogram(output_directory, "feature_multiplicity_histogram.html", "Feature Multiplicity", "Deduplicated UMI Count", "Frequency", feature_hist[0]);
    }
        // After writing the text files, generate the heatmaps
    
    // --- Convert feature_hist (GArray**) to a dense matrix for the heatmap ---
    int max_len = 0;
    for (int i = 1; i <= features->number_of_features; i++) {
        if (feature_hist[i] && feature_hist[i]->len > max_len) {
            max_len = feature_hist[i]->len;
        }
    }
    
    int max_deduped_count = (max_len > 0) ? max_len - 1 : 0;
    
    int **deduped_histograms =
        malloc((features->number_of_features + 1) * sizeof(int *));
    for (int i = 0; i <= features->number_of_features; ++i) {
        deduped_histograms[i] =
            calloc(max_deduped_count + 1, sizeof(int));
    }

    for (int i = 1; i <= features->number_of_features; ++i) {
        if (feature_hist[i]) {
            for (guint j = 0; j < feature_hist[i]->len; j++) {
                deduped_histograms[i][j] = g_array_index(feature_hist[i], uint32_t, j);
            }
        }
    }

    // Now, use the populated deduped_histograms to get total counts and check the condition.
    int num_features_for_deduped_heatmap = 0;
    int temp_total_deduped_counts[features->number_of_features];
    for (int i = 0; i < features->number_of_features; ++i) {
        long feature_total = 0;
        // The histogram at index j stores the frequency of count j.
        // Sum j * frequency(j) to get the total counts for the feature.
        for (int j = 1; j <= max_deduped_count; j++) {
            feature_total += (long)j * deduped_histograms[i+1][j];
        }
        temp_total_deduped_counts[i] = feature_total;

        if (feature_total > min_heatmap) {
            num_features_for_deduped_heatmap++;
        }
    }

    if (num_features_for_deduped_heatmap < 2) {
        fprintf(stdout, "Skipping heatmap generation as there are fewer than 2 features with sufficient counts.\n");
    } else {
        generate_deduped_heatmap(output_directory, features, deduped_histograms, max_deduped_count, temp_total_deduped_counts, min_heatmap);

        // --- Step 2.5: Calculate feature richness histograms for heatmap ---
        // A more direct approach to building the histograms.

        // Step 1: Find max richness and build an inverted index (feature -> barcodes).
        GHashTable *feature_to_barcodes_hash =
                            g_hash_table_new(g_direct_hash, g_direct_equal);
        int global_max_richness = 0;

        GHashTableIter iter;
        gpointer       barcode_key, deduped_hash_value;

        g_hash_table_iter_init(&iter, barcode_to_deduped_hash);
        while (g_hash_table_iter_next(&iter, &barcode_key, &deduped_hash_value)) {
            char barcode[barcode_length + 1];
            code2string((unsigned char *)barcode_key, barcode, barcode_code_length);    
            if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
            if (filtered_barcodes_hash &&
                g_hash_table_lookup(filtered_barcodes_hash, barcode) == NULL)
                continue;

            GHashTable *features_in_barcode = (GHashTable *)deduped_hash_value;
            guint       n_feats_in_bc       = g_hash_table_size(features_in_barcode);

            if (n_feats_in_bc == 0) continue;

            if (n_feats_in_bc > global_max_richness) {
                global_max_richness = n_feats_in_bc;
            }

            GList *feat_idxs = g_hash_table_get_keys(features_in_barcode);
            for (GList *l = feat_idxs; l; l = l->next) {
                int f_idx = GPOINTER_TO_INT(l->data);
                if (f_idx <= 0) continue;

                GList *barcodes_for_feature =
                    g_hash_table_lookup(feature_to_barcodes_hash, l->data);
                barcodes_for_feature = g_list_prepend(barcodes_for_feature, barcode_key);
                g_hash_table_insert(feature_to_barcodes_hash, l->data,
                                    barcodes_for_feature);
            }
            g_list_free(feat_idxs);
        }

        // Step 2: Allocate and zero the dense histogram matrix.
        int **richness_histograms =
            malloc((features->number_of_features + 1) * sizeof(int *));
        for (int i = 0; i <= features->number_of_features; ++i) {
            richness_histograms[i] =
                calloc(global_max_richness + 1, sizeof(int));
        }

        // Step 3: Populate the histogram matrix using the inverted index.
        GList *active_features = g_hash_table_get_keys(feature_to_barcodes_hash);
        for (GList *l = active_features; l; l = l->next) {
            int   f_idx    = GPOINTER_TO_INT(l->data);
            GList *barcodes = g_hash_table_lookup(feature_to_barcodes_hash, l->data);

            for (GList *b_node = barcodes; b_node; b_node = b_node->next) {
                char *      barcode_str = (char *)b_node->data;
                GHashTable *deduped_hash =
                    g_hash_table_lookup(barcode_to_deduped_hash, barcode_str);

                int num_feats = g_hash_table_size(deduped_hash);

                if (num_feats > 0) {
                    richness_histograms[f_idx][num_feats]++;
                }
            }
        }
        g_list_free(active_features);

        // Step 4: Set column 0 to the max richness for each feature.
        for (int i = 1; i <= features->number_of_features; ++i) {
            int row_max = 0;
            for (int j = global_max_richness; j >= 1; --j) {
                if (richness_histograms[i][j] > 0) {
                    row_max = j;
                    break;
                }
            }
            richness_histograms[i][0] = row_max;
        }

        // --- Clean up ---
        // The barcode strings are owned by barcode_to_deduped_hash, so we only
        // need to free the lists themselves.
        GList *f2b_values = g_hash_table_get_values(feature_to_barcodes_hash);
        for (GList *l = f2b_values; l; l = l->next) {
            g_list_free(l->data);
        }
        g_list_free(f2b_values);
        g_hash_table_destroy(feature_to_barcodes_hash);

        generate_heatmap(output_directory, features, richness_histograms);

        for (int i = 0; i <= features->number_of_features; ++i)
            free(richness_histograms[i]);
        free(richness_histograms);
    }

    for (int i = 0; i <= features->number_of_features; ++i) {
        free(deduped_histograms[i]);
    }
    free(deduped_histograms);

    char features_file[LINE_LENGTH];
    sprintf(features_file,"%s/features.txt",output_directory);
    FILE *featuresfp=fopen(features_file,"w");
    if(!featuresfp){fprintf(stderr,"Error opening features file %s\n",features_file);exit(1);}
    for(int idx=0;idx<features->number_of_features;idx++)
        fprintf(featuresfp,"%s\n",features->feature_names[idx]);
    fclose(featuresfp);

    /* --- NEW: open the feature-per-cell CSV --- */
    char fpc_file[LINE_LENGTH];
    sprintf(fpc_file, "%s/feature_per_cell.csv", output_directory);
    FILE *fpcfp = fopen(fpc_file, "w");
    if (fpcfp == NULL) {
        fprintf(stderr, "Error opening feature_per_cell.csv\n");
        exit(1);
    }
    /* header */
    fprintf(fpcfp,
            "barcode,num_features,top_feature_index,total_deduped_umi\n");

    g_hash_table_iter_init(&iter, barcode_to_deduped_hash);
    while (g_hash_table_iter_next(&iter, &barcode_key, &deduped_hash_value)) {
        char barcode[barcode_length + 1];
        code2string((unsigned char *)barcode_key, barcode, barcode_code_length);
        if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);

        if (filtered_barcodes_hash &&
            g_hash_table_lookup(filtered_barcodes_hash, barcode) == NULL){
            skipped_barcodes++;
            continue;
        }
        processed_barcodes++;

        /* ------------------------------------
         *  Data for feature_per_cell.csv
         * ------------------------------------*/
        GHashTable *features_in_barcode = (GHashTable *)deduped_hash_value;
        guint       n_feats_in_bc       = g_hash_table_size(features_in_barcode);

        guint32     total_umi  = 0;
        guint32     top_count  = 0;
        guint32     top_feat   = 0;
        gboolean    tie        = FALSE;

        GHashTableIter feat_iter;
        gpointer      feat_k, feat_v;
        g_hash_table_iter_init(&feat_iter, features_in_barcode);
        while (g_hash_table_iter_next(&feat_iter, &feat_k, &feat_v)) {
            guint32 feat_idx  = GPOINTER_TO_UINT(feat_k);
            guint32 feat_cnt  = GPOINTER_TO_UINT(feat_v);
            total_umi        += feat_cnt;

            if (feat_cnt > top_count) {
                top_count = feat_cnt;
                top_feat  = feat_idx;
                tie       = FALSE;
            } else if (feat_cnt == top_count) {
                tie = TRUE;
            }
        }
        if (tie) top_feat = 0;

        /* write the CSV line */
        fprintf(fpcfp, "%s,%u,%u,%u\n",
                barcode, n_feats_in_bc, top_feat, total_umi);

        /* ------------------------------------ */
    }

    /* --- NEW: close CSV --- */
    fclose(fpcfp);
}

int find_closest_barcodes(unsigned char* code,unsigned char *corrected_codes, unsigned char *indices){
    int number_of_variants=0;
    for (int i=0; i<barcode_length; i++){
        unsigned char corrected_bases[3];
        //check whitelist for variant match
        int nmatches=find_variant_match(code, i, corrected_bases);
        if(nmatches){
            if (nmatches+number_of_variants > max_barcode_mismatches){
                return 0;
            }
            for (int j=0; j<nmatches && number_of_variants+j < max_barcode_mismatches; j++){
                memcpy(corrected_codes+(number_of_variants+j)*barcode_code_length, code, barcode_code_length);
                const int code_index = i / 4;
                corrected_codes[(number_of_variants+j)*barcode_code_length + code_index] = corrected_bases[j];
                indices[number_of_variants+j]=i;
            }
            number_of_variants+=nmatches;
        }
        if (number_of_variants > max_barcode_mismatches){
            return 0;
        }
    }
    return number_of_variants;
}
int find_variant_match(unsigned char *code, int sequence_index, unsigned char *corrected_bases ){
    const unsigned char index= sequence_index / 4;
    const unsigned char shift= 6-2*(sequence_index % 4);
    unsigned char mod_code[barcode_code_length];
    memcpy(mod_code, code, barcode_code_length);
    int number_of_variants=0; 
    for (unsigned char i=1; i<4; i++){
        unsigned char mod = code[index] ^ (i << shift );
        mod_code[index]=mod;
        if (g_hash_table_lookup(whitelist_hash, mod_code) != NULL){
            corrected_bases[number_of_variants++]=mod;
        }
    }
    return number_of_variants;
}
void process_pending_barcodes( data_structures *hashes, memory_pool_collection *pools, statistics *stats, double min_posterior){
    unmatched_barcodes_features_block *current_entry_block=stats->unmatched_list.first_entry;
    while (current_entry_block != NULL){
        unsigned char *retcode=find_best_posterior_match(current_entry_block, number_of_features, min_posterior,stats, hashes);
        if (retcode != 0){
            unmatched_barcodes_features current_entry;
            read_unmatched_barcodes_features_block(current_entry_block, &current_entry);
            char umi[umi_length+1];
            code2string(current_entry.umi,umi, umi_code_length);
            update_feature_counts_from_code(retcode, umi, current_entry.feature_index, hashes,pools);
        }
        current_entry_block=current_entry_block->next;
    }
}
unsigned char* find_best_posterior_match (unmatched_barcodes_features_block *entry_block, int number_of_features, double min_posterior, statistics *stats, data_structures *hashes){
    unmatched_barcodes_features entry_struct;
    unmatched_barcodes_features *entry=&entry_struct;
    read_unmatched_barcodes_features_block(entry_block, entry);
    if (!entry->feature_index){
        return 0;
    }
    double priors[max_barcode_mismatches+1];
    double evidence[max_barcode_mismatches+1];
    double total_evidence = 0.0;
    for (int i=0; i<entry->number_of_closest_barcodes; i++){
        //find counts for the barcode
        int total_counts=1;
        unsigned char* barcode =(entry->closest_barcodes) + i*barcode_code_length;
        feature_counts *s=g_hash_table_lookup(hashes->filtered_hash, barcode);
        if (s != NULL){
            total_counts += GPOINTER_TO_UINT(g_hash_table_lookup(s->counts, GUINT_TO_POINTER(0)));
        }
        priors[i]= pow(10,-0.1*(entry->Qscores[i]-33));
        evidence[i]=(double) total_counts * priors[i];
        total_evidence+=evidence[i];
    }
    for (int i = 0; i < entry->number_of_closest_barcodes; i++) {
        double posteriors = evidence[i] / total_evidence;
        if(posteriors > min_posterior){
            stats->pending_recovered++;
            return (entry->closest_barcodes) + i*barcode_code_length;
        }   
    }
    return 0;
}
int simpleCorrectFeature(char *line, feature_arrays *features, int maxN, int maxHammingDistance, int *hamming_distance){
    int best_feature=simple_search(features, line);
    *hamming_distance=0;
    if (!best_feature){
        const size_t length=strlen(line)-1;
        char buffer[(length+1) * (4 << ((maxN-1)*2))];
        char *corrected_seqs[ 4 << ((maxN-1)*2)];
        int nAlts=checkSequenceAndCorrectForN(line, corrected_seqs, buffer, length, maxN);
        if (nAlts == 1){
            return simple_hamming_search(features, corrected_seqs[0], maxHammingDistance,hamming_distance);
        }
        if (nAlts > 1){
            int bestHammingDistance=maxHammingDistance;
            int ambiguous=0;
            for (int i=0; i<nAlts; i++){
                int feature_index=simple_search(features, corrected_seqs[i]);
                if (feature_index){
                    *hamming_distance=0;
                    return feature_index;
                }
            }
            for (int i=0; i<nAlts; i++){
                int myHammingDistance;
                int feature_index=simple_hamming_search(features, corrected_seqs[i], maxHammingDistance,&myHammingDistance);
                if (feature_index && myHammingDistance < bestHammingDistance){
                    bestHammingDistance=myHammingDistance;
                    best_feature=feature_index;
                    ambiguous=0;
                }
                else if (feature_index && myHammingDistance == bestHammingDistance){
                    if(best_feature){
                        ambiguous=1;
                    }
                    else{
                        best_feature=feature_index;
                        bestHammingDistance=myHammingDistance;
                        ambiguous=0;
                    }
                }
            }
            if(best_feature && !ambiguous){
                *hamming_distance=bestHammingDistance;
                return best_feature;
            }
        }

    }
    return best_feature;
}
int checkAndCorrectFeature(char *line, feature_arrays *features,int maxHammingDistance, int nThreads, int *hamming_distance, char *matching_sequence, int maxN,char *ambiguous, uint16_t *match_position){
    const size_t length=strlen(line)-1;
    char buffer[(length+1) * (4 << ((maxN-1)*2))];
    char *corrected_seqs[ 4 << ((maxN-1)*2)];
    int hamming=0;
    //return ambiguous if the hamming distance is non-zero but the feature is zero
    //return ambiguous if there are too many Ns ie. nAlts is zero - distinguish this by setting hamming distance to maxHammingDistance+1
    int nAlts=checkSequenceAndCorrectForN(line, corrected_seqs, buffer, length, maxN);
    if (!nAlts){
        *hamming_distance=maxHammingDistance+1;
        *ambiguous=1;
        return 0;
    }
    if (nAlts == 1){
        char *myMatchingSequence;
        uint16_t myMatchPosition;
        //no Ns and the barcode is good in terms of ACGT 
        int feature_index=find_feature_match_parallel(features, line, maxHammingDistance,nThreads,&hamming,&myMatchingSequence, &myMatchPosition);
        *hamming_distance=hamming;
        if (feature_index){
            memcpy(matching_sequence, myMatchingSequence, features->feature_lengths[feature_index-1]);
            matching_sequence[features->feature_lengths[feature_index-1]]='\0';
            *match_position = myMatchPosition;
        }
        else if (hamming <= maxHammingDistance){
            *ambiguous=1;
        }    
        return feature_index;
    }
    else{
        //if there are Ns in the barcode then we need to check all the possible sequences
        //and return if there is a unique match in the whitelist
        int myAmbiguous=0;
        int bestfeature_index=0;
        int bestHammingDistance=maxHammingDistance;
        char *bestMatchingSequence=0; 
        uint16_t bestMatchPosition=0;
        //int bestAlt=0;

        for (int i=0; i<nAlts; i++){
            //DEBUG_PRINT( "Checking %s\n", corrected_seqs[i]);
            char *myMatchingSequence;
            uint16_t myMatchPosition;
            int feature_index=find_feature_match_parallel(features, corrected_seqs[i], maxHammingDistance,nThreads,&hamming,&myMatchingSequence, &myMatchPosition);
            //DEBUG_PRINT( "Feature index %d %d\n", feature_index,hamming);
            
            if (feature_index && !hamming){
                memcpy(matching_sequence, myMatchingSequence, features->feature_lengths[feature_index-1]);
                matching_sequence[features->feature_lengths[feature_index-1]]='\0';
                *hamming_distance=hamming;
                *match_position = myMatchPosition;
                return feature_index;
            }
            else if (hamming < bestHammingDistance){
                    bestHammingDistance=hamming;
                    bestfeature_index=feature_index;
                    bestMatchingSequence=myMatchingSequence;
                    myAmbiguous=0;
                    bestMatchPosition=myMatchPosition;
            }
            else if (hamming == bestHammingDistance){
                if( bestfeature_index){
                    myAmbiguous=1;
                }
            }
        }
        if (bestHammingDistance <= maxHammingDistance){
            *hamming_distance=bestHammingDistance;
            *ambiguous=myAmbiguous;
            if (!myAmbiguous && bestfeature_index){
                memcpy(matching_sequence, bestMatchingSequence, features->feature_lengths[bestfeature_index-1]);
                matching_sequence[features->feature_lengths[bestfeature_index-1]]='\0';
                *match_position = bestMatchPosition;
                return bestfeature_index;
            }
            return 0;
        }
        return 0;
    }
}
size_t barcode_code2number(unsigned char *code){
    unsigned char *key=(unsigned char*)g_hash_table_lookup(whitelist_hash, code);
    if (key == NULL){
        return 0;
    }
    return (key-whitelist)/barcode_code_length+1;
}
int checkAndCorrectBarcode(char **lines, int maxN, uint32_t feature_index, uint16_t match_position, data_structures *hashes, memory_pool_collection *pools, statistics *stats, int barcode_constant_offset){
    char *sequence=lines[0]+barcode_constant_offset;
    //DEBUG_PRINT( "Checking barcode %s\n", lines[1]);
    if (strlen(sequence) < barcode_length + umi_length){
        fprintf(stderr, "Error: Incomplete barcode %s\n", lines[0]);
        return 0;
    }
    if (!check_sequence(sequence + barcode_length, umi_length)){
        stats->nMismatches++;
        return 0;
    }
    //The return code should indicate whether we should calculate the barcode or not
    char buffer[(barcode_length + 1) * (4 << ((max_barcode_n-1)*2))];
    char *corrected_seqs[ 4 << ((max_barcode_n-1)*2)];
    char *candidateBarcode = sequence;
    memset (buffer, 0, (barcode_length + 1) * (4 << ((max_barcode_n-1)*2)));
 
    //DEBUG_PRINT( "Checking barcode %s\n", candidateBarcode);    
    int nAlts=checkSequenceAndCorrectForN(candidateBarcode, corrected_seqs, buffer, barcode_length, maxN);
    if (!nAlts){
        stats->nMismatches++;
        return 0;
    }
    unsigned char code[barcode_code_length];
    if (nAlts == 1){
        //no Ns and the barcode is good in terms of ACGT 
        //check if the barcode is in the filtered whitelist 
        string2code(candidateBarcode, barcode_length, code); 
        if (g_hash_table_lookup(whitelist_hash, code)){
            update_feature_counts_from_code(code, sequence+barcode_length, feature_index, hashes, pools);
            stats->valid++;            
            return 1;
        }
        //if the barcode is not in the whitelist then try to correct it 
        unsigned char indices[max_barcode_mismatches+1];
        unsigned char corrected_codes[(max_barcode_mismatches+1)*barcode_code_length];
        stats->nMismatches++;
        //try to correct it by finding closest barcodes
        int nMatches=find_closest_barcodes(code,corrected_codes, indices);
        if (!nMatches){
            return 0;
        }

        //if unique closest barcode found then correct
        if (nMatches==1){
            char corrected_barcode[barcode_length+1];   
            code2string(corrected_codes, corrected_barcode, barcode_code_length);
            memcpy(sequence, corrected_barcode, barcode_length);
            stats->recovered++;
            update_feature_counts_from_code(corrected_codes, sequence+barcode_length, feature_index, hashes, pools);
            stats->valid++;
            return 1;
        }
        else{
            //otherwise store it for later processing when we have counts to calculate the posterior probabilities
            if (nMatches > 1){
                unsigned char qscores[max_barcode_mismatches+1];
                for (int i=0; i<nMatches; i++){
                    qscores[i]=lines[1][indices[i]];
                }
                add_unmatched_barcode_store_feature(code, corrected_codes, sequence+barcode_length, qscores,feature_index,nMatches, match_position, pools,stats);
                stats->pending++;
            }
            return 0;
        }
        return 0;
    }
    else{
        stats->nMismatches++;
        //if there are Ns in the barcode then we need to check all the possible sequences
        //and return if there is a unique match in the whitelist
        int number_of_matches=0;
        unsigned char corrected_code[barcode_code_length];
        for (int i=0; i<nAlts; i++){
            unsigned char code[barcode_code_length];
            string2code(corrected_seqs[i], barcode_length, code);
            char corrected_barcode[barcode_length+1];   
            code2string(code, corrected_barcode, barcode_code_length);
            if (g_hash_table_lookup(whitelist_hash, code)){
                if (number_of_matches){
                    return 0;
                }
                else{
                    number_of_matches++;
                    memcpy(corrected_code, code, barcode_code_length);  
                }
            }
        }
        if (number_of_matches == 1){
            stats->recovered++;
            char corrected_barcode[barcode_length+1];   
            code2string(corrected_code, corrected_barcode, barcode_code_length);
            memcpy(sequence, corrected_barcode, barcode_length);
            update_feature_counts_from_code(corrected_code, sequence + barcode_length,feature_index, hashes,pools);
            stats->valid++;
            return 1;
        }
        return 0;
    }
}


void finalize_processing(feature_arrays *features, data_structures *hashes,  char *directory, memory_pool_collection *pools, statistics *stats, uint16_t stringency, uint16_t min_counts, double min_posterior, GHashTable *filtered_barcodes_hash){
    process_pending_barcodes(hashes, pools, stats,min_posterior);
    double elapsed_time = get_time_in_seconds() - stats->start_time;
    fprintf(stderr, "Finished processing %ld reads in %.2f seconds (%.1f thousand reads/second)\n", stats->number_of_reads, elapsed_time, stats->number_of_reads / (double)elapsed_time / 1000.0);

    //find size of whitelist_hash and non_whitelist_hash
    int total_feature_counts[features->number_of_features];
    int total_deduped_counts[features->number_of_features];
    int total_barcoded_counts[features->number_of_features];

    GArray **feature_hist = g_new0(GArray*, features->number_of_features + 1);

    
    // Create and populate the deduped counts hash table once.
    GHashTable* barcode_to_deduped_hash = g_hash_table_new_full(hash_int32, equal_int32, NULL, (GDestroyNotify)g_hash_table_destroy);
    find_deduped_counts(hashes, barcode_to_deduped_hash, stringency, min_counts);

    // Print the unfiltered results
    printFeatureCounts(features, total_deduped_counts, total_barcoded_counts, feature_hist, directory, hashes, stats, barcode_to_deduped_hash, NULL);
    
    // If a filter is provided, print the filtered results as well.
    if (filtered_barcodes_hash) {
        // Reset co-expression and histogram data before the second ru
        for (int i = 0; i < features->number_of_features + 1; i++) {
            if (feature_hist[i]) {
                g_array_unref(feature_hist[i]);
                feature_hist[i] = NULL;
            }
        }
        
        printFeatureCounts(features, total_deduped_counts, total_barcoded_counts, feature_hist, directory, hashes, stats, barcode_to_deduped_hash, filtered_barcodes_hash);
    }
    
    // Clean up the hash table now that all printing is done.
    g_hash_table_destroy(barcode_to_deduped_hash);
    
    //DEBUG_PRINT( "Number of reads matched to a feature %ld\n", valid);
    print_feature_sequences(features, total_feature_counts, directory, hashes);
    char deduped_hist_filename[FILENAME_LENGTH];
    sprintf(deduped_hist_filename, "%s/deduped_counts_histograms.txt", directory);
    FILE *deduped_hist_fp = fopen(deduped_hist_filename, "w");
    if (deduped_hist_fp == NULL) {
        perror("Failed to open deduped counts histograms file");
        exit(EXIT_FAILURE);
    }
    fprintf(deduped_hist_fp, "FeatureName\tDedupedCount\tFrequency\n");
    for (int i = 0; i < features->number_of_features; i++) {
        GArray *h = feature_hist[i + 1]; // Histograms are 1-indexed
        if (h && h->len > 0) {
            for (unsigned int j = 0; j < h->len; j++) {
                uint32_t frequency = g_array_index(h, uint32_t, j);
                if (frequency > 0) {
                    fprintf(deduped_hist_fp, "%s\t%d\t%u\n", features->feature_names[i], j, frequency);
                }
            }
        }
    }
    fclose(deduped_hist_fp);

    // EM-related output files removed - no longer needed without EM functionality

    int feature_printed[features->number_of_features];
    memset(feature_printed, 0, features->number_of_features * sizeof(int));

    // EM fitting removed - using simple counts only

    // Process cumulative histogram (simplified - no EM fitting)
    {
        long cumulative_deduped_counts = 0;
        long cumulative_barcoded_counts = 0;
        long cumulative_feature_counts = 0;
        for (int i = 0; i < features->number_of_features; i++) {
            cumulative_deduped_counts += total_deduped_counts[i];
            cumulative_barcoded_counts += total_barcoded_counts[i];
            cumulative_feature_counts += total_feature_counts[i];
        }

        // Output summary to stderr for logging
        fprintf(stderr, "Cumulative counts: deduped=%ld, barcoded=%ld, features=%ld\n",
                cumulative_deduped_counts, cumulative_barcoded_counts, cumulative_feature_counts);
        
        // Generate simple histogram plot without EM fit
        plot_simple_histogram(directory, "umi_counts_histogram", "UMI Counts Distribution", 
                             "UMI Count", "Frequency", feature_hist[0]);
    }

    // Free the feature histogram arrays
    for(int i=0; i < features->number_of_features + 1; ++i)
        if(feature_hist[i]) g_array_unref(feature_hist[i]);
    g_free(feature_hist);
}

void open_fastq_files(const char *barcode_fastq, const char *forward_fastq, const char *reverse_fastq, gzFile *barcode_fastqgz, gzFile *forward_fastqgz, gzFile *reverse_fastqgz) {
    *barcode_fastqgz = gzopen(barcode_fastq, "r");
    if (*barcode_fastqgz == NULL) {
        fprintf(stderr, "Error: Unable to open R1 FASTQ file %s\n", barcode_fastq);
        exit(EXIT_FAILURE);
    }
    if (forward_fastq != NULL) {
        *forward_fastqgz = gzopen(forward_fastq, "r");
        if (*forward_fastqgz == NULL) {
            fprintf(stderr, "Error: Unable to open R2 FASTQ file %s\n", forward_fastq);
            exit(EXIT_FAILURE);
        }
    }
    if (reverse_fastq != NULL) {
        *reverse_fastqgz = gzopen(reverse_fastq, "r");
        if (*reverse_fastqgz == NULL) {
            fprintf(stderr, "Error: Unable to open R3 FASTQ file %s\n", reverse_fastq);
            exit(EXIT_FAILURE);
        }
    }
}
fastq_reader* allocate_fastq_reader( char **filenames, int nfiles, int filetype, size_t read_size, size_t read_buffer_lines) {
    fastq_reader *reader = malloc(sizeof(fastq_reader));
    int total_filename_length = 0;
    for (int i = 0; i < nfiles; i++) {
        total_filename_length += strlen(filenames[i])+1; 
    }
    reader->filenames = malloc(nfiles * sizeof(char *));
    reader->concatenated_filenames=malloc(total_filename_length);

    char *this_filename = reader->concatenated_filenames;
    for (int i = 0; i < nfiles; i++) {
        fprintf (stderr, "Filename %s\n", filenames[i]);
        //do no increment the last time to avoid overun of the buffe
    }

    for (int i = 0; i < nfiles; i++) {
        strcpy(this_filename, filenames[i]);
        reader->filenames[i] = this_filename;
        this_filename += strlen(filenames[i]) + 1;
    }
    reader->nfiles = nfiles;
    reader->gz_pointer = NULL;            // Initialize gz_pointer as needed
    reader->filetype = filetype;          // Set the filetype
    return reader;
}
fastq_reader_set *  allocate_fastq_reader_set( char **barcode_filenames, char **forward_filenames, char **reverse_filenames, int nfiles, size_t read_size, size_t read_buffer_lines) {
    fastq_reader_set *reader_set=malloc(sizeof(fastq_reader_set));
    fprintf(stderr, "allocating barcode reader\n");
    reader_set->barcode_reader=allocate_fastq_reader(barcode_filenames, nfiles, 1, read_size, read_buffer_lines);
    reader_set->forward_reader=NULL;
    reader_set->reverse_reader=NULL;
    if (forward_filenames != NULL ) {
        fprintf(stderr, "allocating forward reader\n");
        reader_set->forward_reader=allocate_fastq_reader(forward_filenames, nfiles, 2, read_size, read_buffer_lines);
    }
    if (reverse_filenames != NULL) {
        fprintf(stderr, "allocating reverse reader\n");
        reader_set->reverse_reader=allocate_fastq_reader(reverse_filenames, nfiles, 3, read_size, read_buffer_lines);
    }
    reader_set->read_buffer_lines = read_buffer_lines;
    reader_set->buffer_storage    = malloc(read_buffer_lines * (read_size + 1));
    reader_set->buffer            = malloc(read_buffer_lines * sizeof(char*));
    for (size_t i = 0; i < read_buffer_lines; ++i)
        reader_set->buffer[i] = reader_set->buffer_storage + i * (read_size + 1);

    reader_set->produce_index = reader_set->consume_index = reader_set->filled = 0;
    reader_set->done = 0;

    pthread_mutex_init(&reader_set->mutex, NULL);
    pthread_cond_init (&reader_set->can_produce, NULL);
    pthread_cond_init (&reader_set->can_consume, NULL);
    return reader_set;
}

// Producer thread function
void *read_fastqs_by_set(void *arg) {
    fastq_reader_set *set = (fastq_reader_set *)arg;
    const int thread_id = set->thread_id;
    long long number_of_reads=0;
    const int number_of_readers = (set->forward_reader != NULL) + (set->reverse_reader != NULL) + 1;
    fastq_reader *readers[number_of_readers];
    readers[0]  = set->barcode_reader;
    const int lines_per_block   = 2 * number_of_readers;
    const int read_buffer_lines = set->read_buffer_lines;
    if (set->forward_reader != NULL) {
        readers[1] = set->forward_reader;
    }
    if (set->reverse_reader != NULL) {
        readers[number_of_readers - 1] = set->reverse_reader;
    }

    char line1[number_of_readers][LINE_LENGTH];
    char line2[number_of_readers][LINE_LENGTH];
    char dummy_line[LINE_LENGTH];
    int done = 0;
    // Open the files
    int file_index=0;
    for (int j=0; j<number_of_readers; j++){
        fprintf(stderr, "Opening file %s\n", readers[j]->filenames[file_index]);
        readers[j]->gz_pointer = gzopen(readers[j]->filenames[file_index], "rb");
        if (readers[j]->gz_pointer == NULL) {
            fprintf(stderr, "Error: Unable to open file %s\n", readers[j]->filenames[file_index]);
            exit(EXIT_FAILURE);
        }
    }
    while(!done){
        int eof_detected = 0;
        for (int j = 0; j < number_of_readers; j++) {
            if (gzgets(readers[j]->gz_pointer, dummy_line, LINE_LENGTH) == Z_NULL ||
                gzgets(readers[j]->gz_pointer, line1[j], LINE_LENGTH) == Z_NULL ||
                gzgets(readers[j]->gz_pointer, dummy_line, LINE_LENGTH) == Z_NULL ||
                gzgets(readers[j]->gz_pointer, line2[j], LINE_LENGTH) == Z_NULL)
            {
                if (!gzeof(readers[j]->gz_pointer)) {
                    fprintf(stderr, "Thread %d: Error or incomplete record in FASTQ file %s\n",
                            thread_id, readers[j]->filenames[file_index]);
                    exit(EXIT_FAILURE);
                }
                eof_detected = 1;
                break; // Exit inner loop immediately
            }
        }

        if (eof_detected || (max_reads && number_of_reads >= max_reads)) {
            for (int j = 0; j < number_of_readers; j++) {
                if (readers[j]->gz_pointer) {
                    gzclose(readers[j]->gz_pointer);
                    readers[j]->gz_pointer = NULL;
                }
            }
            if (++file_index >= readers[0]->nfiles) {
                done = 1; // All files for this producer are processed
            } else {
                // Open the next set of files
                for (int j = 0; j < number_of_readers; j++) {
                    fprintf(stderr, "Thread %d opening next file: %s\n", thread_id, readers[j]->filenames[file_index]);
                    readers[j]->gz_pointer = gzopen(readers[j]->filenames[file_index], "rb");
                    if (readers[j]->gz_pointer == NULL) {
                        fprintf(stderr, "Error: Unable to open file %s\n", readers[j]->filenames[file_index]);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            continue; // Restart the while loop
        }

        if (number_of_reads > 0 && number_of_reads % 1000000 == 0) {
            fprintf(stderr, "Thread %d has read %lld million records.\n", thread_id, number_of_reads / 1000000);
        }
        number_of_reads++;

        /* ---------- NEW: write a whole block into set->buffer ---------- */
        pthread_mutex_lock(&set->mutex);

        while (set->filled >= read_buffer_lines - lines_per_block)
            pthread_cond_wait(&set->can_produce, &set->mutex);

        size_t p = set->produce_index;
        /* barcode (always) */
        strcpy(set->buffer[p],                             line1[0]);
        strcpy(set->buffer[(p+1)%read_buffer_lines],       line2[0]);

        int offset = 2;   /* start of next reader's lines */
        int r_idx  = 1;   /* reader index we are copying from */

        if (set->forward_reader) {            /* forward read present */
            strcpy(set->buffer[(p+offset)%read_buffer_lines],     line1[r_idx]);
            strcpy(set->buffer[(p+offset+1)%read_buffer_lines],   line2[r_idx]);
            offset += 2;
            r_idx  += 1;
        }
        if (set->reverse_reader) {            /* reverse read present */
            strcpy(set->buffer[(p+offset)%read_buffer_lines],     line1[r_idx]);
            strcpy(set->buffer[(p+offset+1)%read_buffer_lines],   line2[r_idx]);
        }

        set->produce_index = (p + lines_per_block) % read_buffer_lines;
        set->filled       += lines_per_block;
        pthread_cond_signal(&set->can_consume);
        pthread_mutex_unlock(&set->mutex);
    }
    
    // Signal that this producer is completely finished
    pthread_mutex_lock(&set->mutex);
    set->done = 1;
    pthread_cond_broadcast(&set->can_consume);
    pthread_mutex_unlock(&set->mutex);
    
    fprintf(stderr, "Thread %d done reading\n", thread_id);
    pthread_exit(NULL);
}
void process_multiple_feature_sequences(int nsequences, char **sequences, int *orientations, feature_arrays *features, int maxHammingDistance, int nThreads, int feature_constant_offset, int max_feature_n, uint32_t *feature_index, int *hamming_distance, char *matching_sequence, uint16_t *match_position) {
    //check if the sequences are the same

    int best_hamming_distance = maxHammingDistance;
    uint32_t best_feature_index = 0;
    char best_matching_sequence[LINE_LENGTH];
    uint16_t best_match_position=0;
    char ambiguous = 0;
    int myHammingDistance=maxHammingDistance;
    for (int i=0; i<nsequences; i++){
        DEBUG_PRINT( "Sequence %s\n", sequences[i]);
        if(sequences[i]){
            uint32_t myFeatureIndex=0;
            char my_matching_sequence[LINE_LENGTH];
            uint16_t my_match_position=0;
            if (orientations[i] != 1){
                process_feature_sequence(sequences[i], features, maxHammingDistance, nThreads, feature_constant_offset, max_feature_n, &myFeatureIndex, &myHammingDistance, my_matching_sequence, &my_match_position);
                if (myFeatureIndex && myHammingDistance == best_hamming_distance) {
                    if (best_feature_index) {
                        ambiguous = 1;
                    }
                    else {
                        //first time we have a match
                        ambiguous = 0;
                        best_hamming_distance = myHammingDistance;
                        best_feature_index = myFeatureIndex;
                        best_match_position=my_match_position;
                        strcpy(best_matching_sequence, my_matching_sequence);
                    }
                }
                else if (myFeatureIndex && myHammingDistance < best_hamming_distance) {
                    ambiguous = 0;
                    best_hamming_distance = myHammingDistance;
                    best_feature_index = myFeatureIndex;
                    best_match_position=my_match_position;
                    strcpy(best_matching_sequence, my_matching_sequence);
                }
            }
            if (orientations[i] != 0){
                char reverse_sequence[LINE_LENGTH];
                strcpy(reverse_sequence, sequences[i]);
                reverse_complement_in_place(reverse_sequence);
                process_feature_sequence(reverse_sequence, features, maxHammingDistance, nThreads, feature_constant_offset, max_feature_n, &myFeatureIndex, &myHammingDistance, my_matching_sequence, &my_match_position);
            }
        }
    }
    if (ambiguous || myHammingDistance > maxHammingDistance) {
        *feature_index = 0;
        *hamming_distance = best_hamming_distance;
    } else {
        *feature_index = best_feature_index;
        *hamming_distance = best_hamming_distance;
        *match_position = best_match_position;
        if (best_feature_index) {
            strcpy(matching_sequence, best_matching_sequence);
        }
    }
}
void process_feature_sequence(char *sequence, feature_arrays *features, int maxHammingDistance, int nThreads, int feature_constant_offset, int max_feature_n, uint32_t *feature_index, int *hamming_distance, char *matching_sequence, uint16_t *match_position) {
    if (limit_search != -1 && feature_constant_offset > 0) {
        int bestHammingDistance = maxHammingDistance + 1;
        uint32_t bestFeatureIndex = 0;
        uint16_t bestMatchPosition = 0;
        char ambiguous = 0;

        for (int i = 0; i <= limit_search; i++) {
            // For i=0, this checks offset 0. For i>0, it checks +i and -i.
            for (int sign = 1; sign >= -1; sign -= 2) {
                if (i == 0 && sign == -1) continue; // Only check offset 0 once

                int offset_delta = i * sign;
                int current_offset = feature_constant_offset + offset_delta;
                if (current_offset < 0) continue;

                int currentHammingDistance = 0;
                uint32_t currentFeatureIndex = simpleCorrectFeature(sequence + current_offset, features, max_feature_n, maxHammingDistance, &currentHammingDistance);

                if (currentFeatureIndex) {
                    // If a perfect match is found, we can exit immediately.
                    if (currentHammingDistance == 0) {
                        *feature_index = currentFeatureIndex;
                        *hamming_distance = 0;
                        *match_position = current_offset;
                        memcpy(matching_sequence, sequence + current_offset, features->feature_lengths[currentFeatureIndex - 1]);
                        matching_sequence[features->feature_lengths[currentFeatureIndex - 1]] = '\0';
                        return;
                    }

                    if (currentHammingDistance < bestHammingDistance) {
                        bestHammingDistance = currentHammingDistance;
                        bestFeatureIndex = currentFeatureIndex;
                        bestMatchPosition = current_offset;
                        ambiguous = 0;
                    } else if (currentHammingDistance == bestHammingDistance) {
                        if (bestFeatureIndex != currentFeatureIndex) {
                            ambiguous = 1;
                        }
                    }
                }
            }
        }

        if (!ambiguous && bestFeatureIndex) {
            *feature_index = bestFeatureIndex;
            *hamming_distance = bestHammingDistance;
            *match_position = bestMatchPosition;
            memcpy(matching_sequence, sequence + bestMatchPosition, features->feature_lengths[bestFeatureIndex - 1]);
            matching_sequence[features->feature_lengths[bestFeatureIndex - 1]] = '\0';
        } else {
            *feature_index = 0;
            *hamming_distance = bestHammingDistance;
        }
        return;
    }
    int bestHammingDistance=0;
    int variableMaxHammingDistance=maxHammingDistance;
    uint32_t myFeatureIndex=0;
    if (feature_constant_offset > 0) {
        int constantHammingDistance=0;
        myFeatureIndex = simpleCorrectFeature(sequence + feature_constant_offset, features, max_feature_n, maxHammingDistance, &constantHammingDistance);
        if (myFeatureIndex ) {
            // need to have feature_index -1 because the feature index is 1 based but array in struct is 0 based
            memcpy(matching_sequence, sequence + feature_constant_offset, features->feature_lengths[myFeatureIndex - 1]);
            matching_sequence[features->feature_lengths[myFeatureIndex - 1]] = '\0';
            *match_position = feature_constant_offset;
            bestHammingDistance=constantHammingDistance;
            variableMaxHammingDistance=constantHammingDistance;
            //return immediately if we have a match with distance 0 (no need to check the variable hamming distance because it is the best if it also matches the offset)
            if (constantHammingDistance == 0) {
                *feature_index = myFeatureIndex;
                *hamming_distance = bestHammingDistance;
                return;
            }
        }
        //repeat with feature_constant_offset - 1 and +1 tried higher changes but seem worse
        int minusOne_feature_constant_offset = feature_constant_offset - 1;
        int minusOneFeatureIndex = simpleCorrectFeature(sequence + minusOne_feature_constant_offset, features, max_feature_n, maxHammingDistance, &constantHammingDistance);
        if (minusOneFeatureIndex && constantHammingDistance == 0) {
            *feature_index = minusOneFeatureIndex;
            *hamming_distance = bestHammingDistance;
            *match_position = minusOne_feature_constant_offset;
            return;
        }
        if (minusOneFeatureIndex && constantHammingDistance < bestHammingDistance) {
            bestHammingDistance=constantHammingDistance;
            myFeatureIndex=minusOneFeatureIndex;
            *match_position = minusOne_feature_constant_offset;
        }
        int plusOne_feature_constant_offset = feature_constant_offset + 1;
        int plusOneFeatureIndex = simpleCorrectFeature(sequence + plusOne_feature_constant_offset, features, max_feature_n, maxHammingDistance, &constantHammingDistance);
        if (plusOneFeatureIndex && constantHammingDistance == 0) {
            *feature_index = plusOneFeatureIndex;
            *hamming_distance = bestHammingDistance;
            *match_position = plusOne_feature_constant_offset;
            return;
        }
        if (plusOneFeatureIndex && constantHammingDistance < bestHammingDistance) {
            bestHammingDistance=constantHammingDistance;
            myFeatureIndex=plusOneFeatureIndex;
            *match_position = plusOne_feature_constant_offset;
        }

    }
    if (myFeatureIndex && bestHammingDistance) {
        char ambiguous = 0;
        char new_matching_sequence[LINE_LENGTH];
        uint16_t new_match_position;
        int variableHammingDistance=variableMaxHammingDistance;
        int new_feature_index = checkAndCorrectFeature(sequence, features, variableMaxHammingDistance, nThreads, &variableHammingDistance, new_matching_sequence, max_feature_n, &ambiguous, &new_match_position);
        if (variableHammingDistance < variableMaxHammingDistance) {
            myFeatureIndex = new_feature_index; // even if zero keep the new feature index because that would mean ambiguity
            strcpy(matching_sequence, new_matching_sequence);
            bestHammingDistance=variableHammingDistance;
            *match_position=new_match_position;
        }
    } else if (!myFeatureIndex) {
        char ambiguous = 0;
        char new_matching_sequence[LINE_LENGTH];
        uint16_t new_match_position;
        int variableHammingDistance=maxHammingDistance;
        int newFeatureIndex = checkAndCorrectFeature(sequence, features, variableMaxHammingDistance, nThreads, &variableHammingDistance, new_matching_sequence, max_feature_n, &ambiguous, &new_match_position);
        if (!ambiguous && newFeatureIndex && variableHammingDistance < maxHammingDistance) {
            bestHammingDistance=variableHammingDistance;
            strcpy(matching_sequence, new_matching_sequence);
            myFeatureIndex = newFeatureIndex;
            *match_position=new_match_position;
        }
    }
    if (bestHammingDistance > maxHammingDistance) {
        myFeatureIndex = 0;
    }
    *feature_index = myFeatureIndex;
    *hamming_distance = bestHammingDistance;
}
void *consume_reads(void *arg) {
    fastq_processor *processor_args = (fastq_processor *)arg;
    const int nsets = processor_args->nsets;
    fastq_reader_set **reader_sets = processor_args->reader_sets;
    const sample_args *sample_args = processor_args->sample_args;

    // Use thread-local statistics and hashes
    statistics *stats;  
    data_structures *hashes;  
    memory_pool_collection *pools;  
    const int thread_id = processor_args->thread_id;  
    // default to global (legacy) structures  
    stats  = &sample_args->stats[thread_id];  
    hashes = &sample_args->hashes[thread_id];  
    pools  = sample_args->pools[thread_id];

    int done=0;
    //first one is always barcode
    char *lines_buffer=malloc(6*LINE_LENGTH);
    char *lines[6];
    char done_flags[nsets];
    memset(done_flags, 0, nsets*sizeof(char));
    const int feature_constant_offset = sample_args->feature_constant_offset;
    const int barcode_constant_offset = sample_args->barcode_constant_offset;
    feature_arrays *features = sample_args->features;
    const int maxHammingDistance = sample_args->maxHammingDistance;
    const int nThreads = sample_args->nThreads;
    int nreaders=(reader_sets[0]->forward_reader && reader_sets[0]->reverse_reader)?3:2;
    const int lines_per_block=2*nreaders;           /* NEW – 4 or 6 lines */
    for (int i=0; i<6; i++){
        lines[i]=lines_buffer+i*LINE_LENGTH;
    }
    char **barcode_lines = lines;
    char **forward_lines = 0;
    char **reverse_lines = 0;
    if (nreaders == 3) {
        forward_lines= lines + 2;
        reverse_lines = lines + 4;
    }
    if (nreaders == 2) {
        if (reader_sets[0]->forward_reader) {
            forward_lines = lines + 2;
        }
        else {
            reverse_lines = lines + 2;
        }
    }
    while (!done) {
        // will store 0 type (barcode)
        int data_available = 0;
        for (int i = 0; i < nsets; i++) {
            if (done_flags[i]){
                continue;
            }
            fastq_reader_set *set = reader_sets[i];
            //check the mutexes
            pthread_mutex_lock(&set->mutex);

            while (set->filled < lines_per_block && !set->done)
                pthread_cond_wait(&set->can_consume, &set->mutex);

            if (set->done && set->filled == 0) {
                pthread_mutex_unlock(&set->mutex);
                done_flags[i] = 1;
                continue;
            }

            /* If we fall through, there's either data or it's the end */
            if (set->filled == 0) { // This means set->done must be true
                pthread_mutex_unlock(&set->mutex);
                done_flags[i] = 1;
                continue;
            }
            
            /* we have a data block – copy it out */
            data_available = 1; // We found work
            size_t c = set->consume_index;
            strcpy(barcode_lines[0], set->buffer[c]);
            strcpy(barcode_lines[1], set->buffer[(c+1) % set->read_buffer_lines]);

            if (forward_lines) {
                strcpy(forward_lines[0], set->buffer[(c+2) % set->read_buffer_lines]);
                strcpy(forward_lines[1], set->buffer[(c+3) % set->read_buffer_lines]);
            }
            if (reverse_lines) {
                int off = (nreaders == 3) ? 4 : 2;
                strcpy(reverse_lines[0], set->buffer[(c+off) % set->read_buffer_lines]);
                strcpy(reverse_lines[1], set->buffer[(c+off+1) % set->read_buffer_lines]);
            }
            //signal that the data has been consumed
            set->consume_index = (set->consume_index + lines_per_block) % set->read_buffer_lines;
            set->filled      -= lines_per_block;
            pthread_cond_signal(&set->can_produce);
            pthread_mutex_unlock(&set->mutex);
            //process the data
            char matching_sequence[LINE_LENGTH];
            int hamming_distance = 0;
            uint32_t feature_index = 0;
            uint16_t match_position = 0;
            int missing_flag=0;
            if (forward_lines && reverse_lines) {
                char *sequences[2]={0,0};
                int orientations[2]={0,0};
                sequences[0] = forward_lines[0];
                sequences[1] = reverse_lines[0];
                orientations[0] = 1;
                orientations[1] = 0;
                process_multiple_feature_sequences(2, sequences, orientations, features, maxHammingDistance, nThreads, feature_constant_offset, max_feature_n, &feature_index, &hamming_distance, matching_sequence, &match_position);
            }
            else {
                char *sequence=(forward_lines)?forward_lines[0]:reverse_lines[0];
                process_feature_sequence(sequence, features, maxHammingDistance, nThreads, feature_constant_offset, max_feature_n, &feature_index, &hamming_distance, matching_sequence, &match_position);
            }

            // The mutex is removed from here
            if (feature_index){
                // need to have feature_index -1 because the feature index is 1 based but array in struct is 0 based
                matching_sequence[features->feature_lengths[feature_index - 1]] = '\0';
                insert_feature_sequence(matching_sequence, feature_index, hamming_distance, match_position, hashes, pools);
                checkAndCorrectBarcode(barcode_lines, max_barcode_n, feature_index, match_position, hashes, pools, stats, barcode_constant_offset);
            }
            else{
                missing_flag=1;
            }
            stats->number_of_reads++;
            if (stats->number_of_reads % 1000000 == 0)
            {
                double elapsed_time = get_time_in_seconds() - stats->start_time;
                fprintf(stderr, "Processed %ld million reads in %.1f seconds\n", stats->number_of_reads / 1000000, elapsed_time);
            }
            if (missing_flag){
                stats->nMismatches++;
                stats->total_unmatched_features++;
            }
            break; // Exit the for loop to process the data
        }
        for (int i=0; i<nsets; i++){
            if (!done_flags[i]){
                done=0;
                break;
            }
            done=1;
        }
        if(!done && !data_available){
            sched_yield();
        }
    }
    //free the lines buffer
    free(lines_buffer);
    pthread_exit(NULL);
}
void free_fastq_reader(fastq_reader *reader) {
    if (reader == NULL) return;
    // The gz_pointer is closed by the reader thread itself.
    free(reader->concatenated_filenames);
    free(reader->filenames);
    free(reader); // Free the reader struct itself here.
}
void free_fastq_reader_set(fastq_reader_set *reader_set) {
    free_fastq_reader(reader_set->barcode_reader);
    if (reader_set->forward_reader) {
        free_fastq_reader(reader_set->forward_reader);
    }
    if (reader_set->reverse_reader) {
        free_fastq_reader(reader_set->reverse_reader);
    }
    free(reader_set->buffer);
    free(reader_set->buffer_storage);
    pthread_mutex_destroy(&reader_set->mutex);
    pthread_cond_destroy(&reader_set->can_produce);
    pthread_cond_destroy(&reader_set->can_consume);
    free(reader_set);
}
void free_fastq_processor(fastq_processor *processor_args) {
    for (int i = 0; i < processor_args->nsets; i++) {
        if (processor_args->reader_sets[i]) {
            free_fastq_reader_set(processor_args->reader_sets[i]);
        }
    }
    pthread_mutex_destroy(&processor_args->process_mutex);
}
void merge_stats(statistics *merged_stats, statistics *thread_stats) {
    merged_stats->nMismatches += thread_stats->nMismatches;
    merged_stats->recovered += thread_stats->recovered;
    merged_stats->pending += thread_stats->pending;
    merged_stats->valid += thread_stats->valid;
    merged_stats->pending_recovered += thread_stats->pending_recovered;
    merged_stats->total_unmatched_features += thread_stats->total_unmatched_features;
    merged_stats->number_of_reads += thread_stats->number_of_reads;
}

typedef struct {
    GHashTable *dst_hash;
    memory_pool_collection *dst_pool;
} merge_context;

static void merge_uint_counters(gpointer key,gpointer value,gpointer user_data)
{
    GHashTable *dest = (GHashTable *)user_data;
    uintptr_t current = GPOINTER_TO_UINT(g_hash_table_lookup(dest, key));
    uintptr_t add     = GPOINTER_TO_UINT(value);
    g_hash_table_replace(dest, key, GUINT_TO_POINTER(current + add));
}

static void copy_uint_entry(gpointer key, gpointer value, gpointer user_data) {
    GHashTable *dest = (GHashTable *)user_data;
    g_hash_table_insert(dest, key, value);
}

void merge_feature_counts(gpointer key, gpointer value, gpointer user_data)
{
    merge_context *ctx = (merge_context *)user_data;
    GHashTable *dst = ctx->dst_hash;
    feature_counts *src_entry = (feature_counts *)value;
    feature_counts *dst_entry = g_hash_table_lookup(dst, key);

    if (dst_entry) {
        /* merge the per-barcode feature counters */
        g_hash_table_foreach(src_entry->counts,
                             merge_uint_counters,
                             dst_entry->counts);
    } else {
        dst_entry = (feature_counts*) allocate_memory_from_pool(ctx->dst_pool->feature_counts_pool);
        memcpy(dst_entry->sequence_code, src_entry->sequence_code, barcode_code_length);
        dst_entry->counts = g_hash_table_new(g_direct_hash, g_direct_equal);
        g_hash_table_foreach(src_entry->counts, copy_uint_entry, dst_entry->counts);
        g_hash_table_insert(dst, dst_entry->sequence_code, dst_entry);
    }
}

void merge_feature_umi_counts(gpointer key, gpointer value, gpointer user_data)
{
    merge_context *ctx        = (merge_context *)user_data;
    GHashTable *dst           = ctx->dst_hash;
    feature_umi_counts *src_ent = (feature_umi_counts *)value;
    feature_umi_counts *dst_ent = g_hash_table_lookup(dst, key);

    if (dst_ent) {
        /* Same barcode-UMI already present – add the counters */
        g_hash_table_foreach(src_ent->counts,merge_uint_counters, dst_ent->counts);
    } else {
        /* Key not present – copy the whole struct into the dst pool and table  */
        dst_ent = (feature_umi_counts*) allocate_memory_from_pool(ctx->dst_pool->feature_umi_counts_pool);
        memcpy(dst_ent->sequence_umi_code, src_ent->sequence_umi_code, 8);
        dst_ent->counts = g_hash_table_new(g_direct_hash, g_direct_equal);
        g_hash_table_foreach(src_ent->counts, copy_uint_entry, dst_ent->counts);
        g_hash_table_insert(dst, dst_ent->sequence_umi_code, dst_ent);
    }
}

void merge_feature_sequences(gpointer key, gpointer value, gpointer user_data) {
    merge_context *ctx = (merge_context *)user_data;
    GHashTable *merged_hash = ctx->dst_hash;
    feature_sequences *thread_entry = (feature_sequences *)value;
    feature_sequences *merged_entry = g_hash_table_lookup(merged_hash, key);
    if (merged_entry) {
        merged_entry->counts += thread_entry->counts;
    } else {
        merged_entry = (feature_sequences*) allocate_memory_from_pool(ctx->dst_pool->feature_sequences_pool);
        memcpy(merged_entry, thread_entry, dynamic_struct_sizes.feature_sequences);
        g_hash_table_insert(merged_hash, merged_entry->sequence, merged_entry);
    }
}

void merge_queues(Queue *dest_q, Queue *src_q) {
    if (!src_q || !src_q->data) return;
    
    while(!is_empty(src_q)){
        enqueue(dest_q, dequeue(src_q));
    }
}

void merge_unmatched_barcodes(unmatched_barcodes_features_block_list *merged_list, unmatched_barcodes_features_block_list *thread_list, memory_pool_collection *merged_pool) {
    unmatched_barcodes_features_block *current_block = thread_list->first_entry;
    while (current_block) {
        unmatched_barcodes_features_block *new_block = (unmatched_barcodes_features_block*)allocate_memory_from_pool(merged_pool->unmatched_barcodes_features_block_pool);
        
        // Copy the data, excluding the 'next' pointer which will be set manually.
        memcpy(new_block, current_block, dynamic_struct_sizes.unmatched_barcodes_features_block);
        new_block->next = NULL;

        // Append to merged list
        if (merged_list->last_entry) {
            merged_list->last_entry->next = new_block;
            merged_list->last_entry = new_block;
        } else {
            merged_list->first_entry = new_block;
            merged_list->last_entry = new_block;
        }
        current_block = current_block->next;
    }

    // Detach the thread list
    thread_list->first_entry = NULL;
    thread_list->last_entry = NULL;
}
void process_files_in_sample(sample_args *args) {
    //allocate buffers here
    //number of lines to read into the buffer
    double  min_posterior=args->min_posterior;
    // gposterior removed - EM functionality no longer needed
    const int sample_index = args->sample_index;
    const int nconsumers = args->consumer_threads_per_set;
    fastq_processor processor_args[nconsumers]; // Array of processor args

    fastq_files_collection *fastq_files=args->fastq_files;
    const int sample_offset=fastq_files->sample_offsets[sample_index];
    const int sample_size=fastq_files->sample_sizes[sample_index];
    

    if (args->filtered_barcodes_name && !args->filtered_barcodes_hash) {
        fprintf(stderr, "Looking for filtered barcodes file at %s\n", args->filtered_barcodes_name);
        char filtered_path[FILENAME_LENGTH];
        snprintf(filtered_path, FILENAME_LENGTH, "%s/%s", args->directory, args->filtered_barcodes_name);
        if (file_exists(filtered_path)) {
            args->filtered_barcodes_hash = g_hash_table_new(g_str_hash, g_str_equal);
            read_barcodes_into_hash(filtered_path, args->filtered_barcodes_hash);
        }
        else {
            fprintf(stderr, "Error: Filtered barcodes file not found at %s\n", filtered_path);
        }
    }


    // Allocate and initialize arrays of statistics and data_structures
    args->stats = malloc(nconsumers * sizeof(statistics));
    args->hashes = malloc(nconsumers * sizeof(data_structures));
    args->pools = malloc(nconsumers * sizeof(memory_pool_collection*));
    if (!args->stats || !args->hashes || !args->pools) {
        perror("Failed to allocate memory for thread-local data");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < nconsumers; i++) {
        initialize_statistics(&args->stats[i]);
        initialize_data_structures(&args->hashes[i]);
        // NEW: Initialize a memory pool for each consumer thread
        args->pools[i] = initialize_memory_pool_collection();
    }

    // Initialize the data structures
    fastq_reader_set *reader_sets[sample_size];
    for (int i = 0; i < sample_size; ++i) {
        char **barcode_files = fastq_files->barcode_fastq + sample_offset + i;
        char **forward_files = (fastq_files->forward_fastq)
                               ? fastq_files->forward_fastq + sample_offset + i
                               : NULL;
        char **reverse_files = (fastq_files->reverse_fastq)
                               ? fastq_files->reverse_fastq + sample_offset + i
                               : NULL;

        reader_sets[i] = allocate_fastq_reader_set(barcode_files ,
                                                   forward_files ,
                                                   reverse_files ,
                                                   /* nfiles = */ 1,
                                                   args->average_read_length,
                                                   args->read_buffer_lines);
        reader_sets[i]->thread_id = i;
    }

    for (int i = 0; i < nconsumers; ++i) {
        processor_args[i].sample_args  = args;
        processor_args[i].reader_sets  = reader_sets;   /* array you built */
        processor_args[i].nsets        = sample_size;   /* one entry per lane */
        processor_args[i].thread_id    = i;
        processor_args[i].nreaders     = (fastq_files->forward_fastq &&
                                           fastq_files->reverse_fastq) ? 3 : 2;

        pthread_mutex_init(&processor_args[i].process_mutex, NULL);
    }

    pthread_t *producer_threads = malloc(sample_size * sizeof(pthread_t));
    for (int i = 0; i < sample_size; ++i) {
        if (pthread_create(&producer_threads[i],
                           NULL,
                           read_fastqs_by_set,
                           reader_sets[i]) != 0)
        {
            perror("Failed to create producer thread");
            exit(EXIT_FAILURE);
        }
    }
    pthread_t consumer_threads[nconsumers];
    fprintf(stderr, "Will use %d threads to process the reads\n", nconsumers);
    for (int j=0; j<nconsumers; j++){
        if (pthread_create(&consumer_threads[j], NULL, consume_reads, (void *)&processor_args[j]) != 0) {
            perror("Failed to create consumer thread");
            exit(EXIT_FAILURE);
        }
    }
    //join the threads
    for (int i = 0; i < sample_size; ++i)
        pthread_join(producer_threads[i], NULL);
    free(producer_threads);
    for (int j=0; j<nconsumers; j++){
        pthread_join(consumer_threads[j], NULL);
    }
    // Merge data from all threads into the first thread's data structures
    for (int i = 1; i < nconsumers; i++) {
        merge_stats(&args->stats[0], &args->stats[i]);
        
        merge_context ctx;
        ctx.dst_pool = args->pools[0];

        ctx.dst_hash = args->hashes[0].filtered_hash;
        g_hash_table_foreach(args->hashes[i].filtered_hash, merge_feature_counts, &ctx);

        ctx.dst_hash = args->hashes[0].sequence_umi_hash;
        g_hash_table_foreach(args->hashes[i].sequence_umi_hash, merge_feature_umi_counts, &ctx);
        
        ctx.dst_hash = args->hashes[0].unique_features_match;
        g_hash_table_foreach(args->hashes[i].unique_features_match, merge_feature_sequences, &ctx);

        merge_unmatched_barcodes(&args->stats[0].unmatched_list, &args->stats[i].unmatched_list, args->pools[0]);
        merge_queues(args->hashes[0].neighbors_queue, args->hashes[i].neighbors_queue);

        // After deep copying all data from thread i, we can free its data structures
        // and its memory pool, as they are no longer needed.
        destroy_data_structures(&args->hashes[i]);
        free_memory_pool_collection(args->pools[i]);
        //[i] = NULL; // Avoid double-free in later cleanup
    }
    // Since merging is not required, finalize using the first thread's data.
    finalize_processing(args->features, &args->hashes[0], args->directory, args->pools[0], &args->stats[0], args->stringency, args->min_counts, min_posterior, args->filtered_barcodes_hash);
   
    // Free the reader sets
    for (int i = 0; i < sample_size; ++i)
        free_fastq_reader_set(reader_sets[i]);
    // need to write a free function for reader_sets
    for (int i = 0; i < nconsumers; i++) {
        pthread_mutex_destroy(&processor_args[i].process_mutex);
    }
    
    // Now that all processing is complete, clean up the resources from thread 0.
    destroy_data_structures(&args->hashes[0]);
    free_memory_pool_collection(args->pools[0]);
    
    // Now free the arrays themselves
    free(args->stats);
    free(args->hashes);
    free(args->pools);
    args->stats = NULL;
    args->hashes = NULL;
    args->pools = NULL;
}

void initialize_data_structures(data_structures *hashes){
    hashes->filtered_hash = g_hash_table_new_full(hash_int32, equal_int32, NULL, destroy_feature_counts);
    hashes->unique_features_match = g_hash_table_new_full(g_str_hash, g_str_equal, NULL, NULL); // Keys/Values are memory-pooled
    // Use g_hash_table_new_full to provide our custom value destroyer function.
    hashes->sequence_umi_hash = g_hash_table_new_full(hash_int64, equal_int64, NULL, destroy_feature_umi_counts);
    
    hashes->neighbors_queue=malloc(sizeof(Queue));
    //check if the memory allocation was successful
    if (hashes->neighbors_queue == NULL) {
        perror("Failed to allocate memory for neighbors queue");
        exit(EXIT_FAILURE);
    }
    init_queue(hashes->neighbors_queue);
}
void destroy_data_structures(data_structures *hashes){
    if (hashes->neighbors_queue){
        if (hashes->neighbors_queue->data)free_queue(hashes->neighbors_queue);
        free(hashes->neighbors_queue);
    }
    if (hashes->filtered_hash) g_hash_table_destroy(hashes->filtered_hash);
    if (hashes->unique_features_match) g_hash_table_destroy(hashes->unique_features_match);
    if (hashes->sequence_umi_hash) g_hash_table_destroy(hashes->sequence_umi_hash);
}


int existing_output_skip(char keep_existing, char *directory){
    if (keep_existing && file_exists(directory)){
        char matrix_filename[4096];
        strcpy(matrix_filename, directory);
        strcat(matrix_filename, "matrix.mtx");
        if (file_exists(matrix_filename)){
            fprintf(stderr, "Matrix file %s found and skipping %s\n", matrix_filename, get_basename(directory));
            return 1;
        }
    }
    return 0;
}

int compare_file_sizes(const void *a, const void *b, void *context) {
    size_t *sizes = (size_t *)context;
    int index_a = *(int *)a;
    int index_b = *(int *)b;

    // Compare sizes in descending order (so larger sizes come first)
    if (sizes[index_a] < sizes[index_b]) {
        return 1;
    } else if (sizes[index_a] > sizes[index_b]) {
        return -1;
    } else {
        return 0;
    }
}

// Custom qsort that allows passing a context (like the sizes array)
void qsort_with_context(void *base, size_t num, size_t size, 
                        int (*compar)(const void *, const void *, void *), 
                        void *context) {
    char *ptr = (char *)base;
    for (size_t i = 0; i < num - 1; i++) {
        for (size_t j = i + 1; j < num; j++) {
            if (compar(ptr + i * size, ptr + j * size, context) > 0) {
                // Swap the elements
                char temp[size];
                memcpy(temp, ptr + i * size, size);
                memcpy(ptr + i * size, ptr + j * size, size);
                memcpy(ptr + j * size, temp, size);
            }
        }
    }
}

// Main function to sort the indices of samples by their size and store the order in sample_order
void cleanup_sample(memory_pool_collection *pools, data_structures *hashes){
    destroy_data_structures(hashes);
    free_memory_pool_collection(pools);
}

void reverse_in_place(char *str) {
    if (str) {
        char *end = str + strlen(str) - 1;
        while (str < end) {
            char temp = *str;
            *str = *end;
            *end = temp;
            str++;
            end--;
        }
    }
}
char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return base;  // For unexpected characters, return the base itself
    }
}

// Function to calculate the reverse complement of a sequence in place
void reverse_complement_in_place(char *seq) {
    int left = 0;
    int right = strlen(seq) - 1;

    // Loop through the sequence from both ends, swapping and complementing
    while (left <= right) {
        // Get the complement of both ends
        char left_complement = complement(seq[left]);
        char right_complement = complement(seq[right]);

        // Swap the left and right complements
        seq[left] = right_complement;
        seq[right] = left_complement;

        left++;
        right--;
    }
}

void reverse_complement_sequence(char *sequence,  char *reverse, int length){
    for (int i = 0; i < length; i++) {
        const char comp_base =match[(unsigned char)sequence[length - i - 1]];
        if (comp_base == 0){
            fprintf(stderr, "Error: Invalid base %c in sequence %s\n", sequence[length - i - 1], sequence);
            exit(EXIT_FAILURE);
        }
        reverse[i] = comp_base;
    }
    reverse[length] = '\0';
}

static char sample_barcode_buffer[LINE_LENGTH];

char* grab_sample_barcode(char *read, uint16_t feature_end, const sample_args *args){
    if (!args->sample_barcodes) return NULL;
    int offset;
    if (args->sample_constant_offset >= 0) offset = args->sample_constant_offset;
    else if (args->sample_offset_relative != 0) offset = (int)feature_end + args->sample_offset_relative;
    else return NULL;
    const int len = args->sample_barcodes->common_length;
    if (offset < 0) return NULL;
    int read_len = strlen(read);
    if (offset + len > read_len) return NULL;
    memcpy(sample_barcode_buffer, read + offset, len);
    sample_barcode_buffer[len] = '\0';
    return sample_barcode_buffer;
}

int find_sample_index(char *seq, const sample_args *args, int *hamming){
    if (!args->sample_barcodes) return 0;
    if (!seq) return 0;
    int hd=0;
    char dummy_match[LINE_LENGTH];
    uint16_t dummy_pos=0;
    int idx = checkAndCorrectFeature(seq, args->sample_barcodes, args->sample_max_hamming, /*threads*/1, &hd, dummy_match, args->sample_max_N, NULL, &dummy_pos);
    if (hamming) *hamming = hd;
    if (idx==0) return 0; // undetermined
    return idx; // 1-based index into sample_barcodes
}





void translate_nxt_inplace(char *barcode, int len) {
    if (!barcode) return;
    if (len >= 9) {
        barcode[7] = complement(barcode[7]);
        barcode[8] = complement(barcode[8]);
    }
}