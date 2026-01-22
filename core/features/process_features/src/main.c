#include "../include/common.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/utils.h"
#include "../include/io.h"
#include <stdio.h>

static void print_usage(const char *prog){
    fprintf(stderr, "\nUsage: %s [options] <FASTQ directories or files>\n\n", prog);
    fprintf(stderr, "Required flags:\n");
    fprintf(stderr, "  -w, --whitelist <file>            10x-style barcode whitelist (one per line)\n");
    fprintf(stderr, "  -f, --featurelist <file>          CSV with 'name' and 'sequence' columns\n");
    fprintf(stderr, "  -d, --directory  <path>           Output directory; one subdir per sample\n\n");

    fprintf(stderr, "Input & Output Files:\n");
    fprintf(stderr, "      --barcode_fastqs    <list>    Comma-separated R1 FASTQ files\n");
    fprintf(stderr, "      --forward_fastqs    <list>    Comma-separated R2 FASTQ files\n");
    fprintf(stderr, "      --reverse_fastqs    <list>    Comma-separated R3 FASTQ files\n");
    fprintf(stderr, "      --barcode_fastq_pattern <str> Pattern to identify barcode FASTQs (default _R1_)\n");
    fprintf(stderr, "      --forward_fastq_pattern <str> Pattern to identify forward FASTQs (default _R2_)\n");
    fprintf(stderr, "      --reverse_fastq_pattern <str> Pattern to identify reverse FASTQs (default _R3_)\n");
    fprintf(stderr, "      --filtered_barcodes <file>    File with barcodes to process, one per line\n");
    fprintf(stderr, "  -k, --keep_existing               Skip processing if output files exist\n");
    fprintf(stderr, "  -a, --as_named                    Treat all input files as single sample\n\n");

    fprintf(stderr, "Barcode & Feature Processing:\n");
    fprintf(stderr, "  -b, --barcode_length    <int>     Length of cell barcode (default 16)\n");
    fprintf(stderr, "  -u, --umi_length        <int>     Length of UMI (default 12)\n");
    fprintf(stderr, "  -o, --feature_constant_offset <int> Expected start position of feature in read (default 0)\n");
    fprintf(stderr, "  -B, --barcode_constant_offset <int> Start position of barcode and UMI (default 0)\n");
    fprintf(stderr, "      --limit_search      <int>     Limit feature search to N bases around offset (-1 = entire read)\n");
    fprintf(stderr, "  -r, --reverse_complement_whitelist Reverse complement whitelist barcodes\n\n");

    fprintf(stderr, "Error Correction & Thresholds:\n");
    fprintf(stderr, "  -m, --maxHammingDistance <int>    Max Hamming distance for feature match (default 1)\n");
    fprintf(stderr, "  -s, --stringency        <int>     UMI dedup stringency (see README, default 1)\n");
    fprintf(stderr, "  -i, --min_counts        <int>     Min reads in UMI clique for counting (default 1)\n");
    fprintf(stderr, "  -M, --min_posterior     <float>   Min posterior probability for barcode rescue (default 0.975)\n");
    fprintf(stderr, "      --max_barcode_mismatches <int> Max mismatches to rescue sequence barcode (default 3)\n");
    fprintf(stderr, "      --feature_n         <int>     Max 'N' bases allowed in feature sequence (default 3)\n");
    fprintf(stderr, "      --barcode_n         <int>     Max 'N' bases allowed in sequence barcode (default 1)\n");
    fprintf(stderr, "      --max_reads         <long>    Max reads to process per FASTQ (0 = all)\n");
    fprintf(stderr, "      --min_prediction    <int>     Min prediction threshold for feature assignment (default 1)\n");
    fprintf(stderr, "      --min_heatmap       <int>     Min deduped count for feature in heatmap (default 0)\n\n");

    fprintf(stderr, "Performance & Parallelism:\n");
    fprintf(stderr, "  -t, --threads           <int>     Max concurrent sample processes (default 8)\n");
    fprintf(stderr, "  -S, --search_threads    <int>     Threads per consumer for feature search (default 4)\n");
    fprintf(stderr, "  -c, --consumer_threads_per_set <int> Consumer threads per sample (default 1)\n");
    fprintf(stderr, "  -R, --read_buffer_lines <int>     Lines for read buffer (default 1024)\n");
    fprintf(stderr, "  -L, --average_read_length <int>   Estimated avg read length for buffer allocation (default 300)\n\n");

    fprintf(stderr, "Miscellaneous:\n");
    fprintf(stderr, "      --translate_NXT               Complement positions 8 and 9 of cell barcodes at output/filter stages\n");
    fprintf(stderr, "  -v, --debug                       Enable verbose debug output\n");
    fprintf(stderr, "  -h, --help                        Show this help and exit\n\n");
}

int main(int argc, char *argv[])
{   
    omp_set_nested(1);
    initseq2Code();
    initcode2seq();
    initdiff2hamming(diff2Hamming);
    initialize_complement();
    feature_code_hash = kh_init(codeu32);
    feature_arrays *features=0;
    int reverse_complement_whitelist=0;
    char *barcodeFastqFilesString=0;
    char *forwardFastqFilesString=0;
    char *reverseFastqFilesString=0;
    char *filtered_barcodes_filename=0;
    char barcode_pattern[LINE_LENGTH], forward_pattern[LINE_LENGTH], reverse_pattern[LINE_LENGTH];
    strcpy(barcode_pattern, "_R1_");
    strcpy(forward_pattern, "_R2_");
    strcpy(reverse_pattern, "_R3_");
    
    int maxHammingDistance=1;
    char directory[LINE_LENGTH];
    char whitelist_filename[4096];
    char sample_flag=1;
    char keep_existing=0;
    uint16_t stringency=1;
    uint16_t min_counts=1;
    int read_buffer_lines=READ_BUFFER_LINES;
    int average_read_length=AVERAGE_READ_LENGTH;
    int feature_constant_offset=0;
    int barcode_constant_offset=0;
    double min_posterior=MIN_POSTERIOR;

    int max_concurrent_processes=8;
    int consumer_threads_per_set=1;
    int search_threads_per_consumer=4;
    int set_consumer_threads_per_set=0;
    int set_search_threads_per_consumer=0;
    uint16_t min_prediction = 1;

    char *sample_barcodes_filename = NULL;
    int sample_max_hamming = 1;
    int sample_max_N = 0;
    int sample_constant_offset_cli = -2;  /* sentinel */
    int sample_offset_relative_cli = 0;
    feature_arrays *sample_barcodes = NULL;

    static struct option long_options[] = {
        {"whitelist", required_argument, 0, 'w'},
        {"featurelist", required_argument, 0, 'f'},
        {"maxHammingDistance", required_argument, 0, 'm'},
        {"feature_constant_offset", required_argument, 0, 'o'},
        {"barcode_constant_offset", required_argument, 0, 'B'},
        {"threads", required_argument, 0, 't'},
        {"search_threads", required_argument, 0, 'S'},
        {"stringency", required_argument, 0, 's'},
        {"min_counts", required_argument, 0, 'i'},
        {"umi_length",required_argument , 0, 'u'},
        {"barcode_length",required_argument , 0, 'b'},
        {"directory", required_argument, 0, 'd'},
        {"debug", no_argument, 0, 'v'},
        {"as_named", no_argument, 0, 'a'},
        {"reverse_complement_whitelist", no_argument, 0, 'r'},
        {"keep_existing", no_argument, 0, 'k'},
        {"read_buffer_lines", required_argument, 0, 'R'},
        {"average_read_length", required_argument, 0, 'L'},
        {"min_posterior", required_argument, 0, 'M'},
        {"consumer_threads_per_set", required_argument, 0, 'c'},
        {"barcode_fastqs", required_argument, 0, 0},
        {"forward_fastqs", required_argument, 0, 1},
        {"reverse_fastqs", required_argument, 0, 2},
        {"max_barcode_mismatches", required_argument, 0, 3},
        {"feature_n", required_argument, 0, 4},
        {"barcode_n", required_argument, 0, 5},
        {"barcode_fastq_pattern", required_argument, 0, 6},
        {"forward_fastq_pattern", required_argument, 0, 7},
        {"reverse_fastq_pattern", required_argument, 0, 8},
        {"max_reads", required_argument, 0, 9},
        {"limit_search", required_argument, 0, 10},
        {"filtered_barcodes", required_argument, 0, 12},
        {"min_prediction", required_argument, 0, 15},
        {"min_heatmap", required_argument, 0, 16},
        {"translate_NXT", no_argument, 0, 17},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    if (argc == 1) { print_usage(argv[0]); return 0; }

    while ((c = getopt_long(argc, argv, "w:b:f:m:s:S:i:t:T:o:d:u:c:vakrDB:R:L:M:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'w': strcpy(whitelist_filename, optarg); break;
            case 'b': barcode_length=atoi(optarg); barcode_code_length=(barcode_length+3)/4; break;
            case 'f': features=read_features_file(optarg); number_of_features=features->number_of_features; maximum_feature_length=features->max_length; feature_code_length=(maximum_feature_length+3)/4; break;
            case 'm': maxHammingDistance=atoi(optarg); break;
            case 's': stringency=(uint16_t)atoi(optarg); break;
            case 'S': set_search_threads_per_consumer=atoi(optarg); break;
            case 'i': min_counts=(uint16_t)atoi(optarg); break;
            case 'd': strcpy(directory, optarg); if (directory[strlen(directory)-1] != '/'){ strcat(directory, "/"); } break;
            case 'o': feature_constant_offset=atoi(optarg); break;
            case 't': max_concurrent_processes=atoi(optarg); break;
            case 'u': umi_length=atoi(optarg); umi_code_length=(umi_length+3)/4; break;
            case 'c': set_consumer_threads_per_set=atoi(optarg); break;
            case 'v': debug = 1; break;
            case 'a': sample_flag=0; break;
            case 'k': keep_existing=1; break;
            case 'r': reverse_complement_whitelist=1; break;
            case 'B': barcode_constant_offset=atoi(optarg); break;
            case 'R': read_buffer_lines=atoi(optarg); break;
            case 'L': average_read_length=atoi(optarg); break;
            case 'M': min_posterior=(double)atof(optarg); break;
            case 0: barcodeFastqFilesString = strdup(optarg); break;
            case 1: forwardFastqFilesString = strdup(optarg); break;
            case 2: reverseFastqFilesString = strdup(optarg); break;
            case 3: max_barcode_mismatches=atoi(optarg); break;    
            case 4: max_feature_n=atoi(optarg); break;
            case 5: max_barcode_n=atoi(optarg); break;
            case 6: strcpy(barcode_pattern, optarg); break;
            case 7: strcpy(forward_pattern, optarg); break;
            case 8: strcpy(reverse_pattern, optarg); break;
            case 9: max_reads=atoll(optarg); break;
            case 10: limit_search = atoi(optarg); break;
            case 12: filtered_barcodes_filename = strdup(optarg); break;
            case 15: min_prediction = atoi(optarg); break;
            case 16: min_heatmap = atoi(optarg); break;
            case 17: translate_NXT = 1; fprintf(stderr, "translate_NXT enabled: complementing positions 8 and 9 at output/filter time.\n"); break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }
    khash_t(strptr) *filtered_barcodes_hash = NULL;
    if (filtered_barcodes_filename) {
        int filtered_barcodes_found = 0;
        if (!file_exists(filtered_barcodes_filename)) {
            printf("filtered_barcodes_filename: %s does not exist  \n", filtered_barcodes_filename);
            char full_path[2048];
            if (strlen(directory) > 0) {
                snprintf(full_path, sizeof(full_path), "%s%s", directory, filtered_barcodes_filename);
                printf("Will check for filtered barcodes file at %s\n", full_path);
                if (file_exists(full_path)) {
                    free(filtered_barcodes_filename);
                    filtered_barcodes_filename = strdup(full_path);
                    filtered_barcodes_found = 1;
                } else {
                    fprintf(stderr, "Warning: Filtered barcodes file not found at %s or %s\n", filtered_barcodes_filename, full_path);
                }
            } else {
                fprintf(stderr, "Warning: Filtered barcodes file not found at %s\n", filtered_barcodes_filename);
            }
        }
        else{
            filtered_barcodes_found = 1;
            printf("Will use filtered barcodes file: %s\n", filtered_barcodes_filename);
        }
        
        if (filtered_barcodes_found) {
            filtered_barcodes_hash = kh_init(strptr);
            read_barcodes_into_hash(filtered_barcodes_filename, filtered_barcodes_hash);
        }
    }
    int positional_arg_count = argc - optind;
    if (positional_arg_count > 0 && (barcodeFastqFilesString != NULL || forwardFastqFilesString != NULL || reverseFastqFilesString != NULL)) {
        fprintf(stderr, "Error: Cannot specify both positional arguments and --barcode_fastqs\n");
        return 1;
    }
    if (strcmp(barcode_pattern, forward_pattern) == 0 || strcmp(barcode_pattern, reverse_pattern) == 0 || strcmp(forward_pattern, reverse_pattern) == 0){
        fprintf(stderr, "Error: Barcode, forward, and reverse patterns must be different\n");
        return 1;
    }
    


    fastq_files_collection fastq_files;
    memset(&fastq_files, 0, sizeof(fastq_files));
    if (positional_arg_count && is_directory(argv[optind])){
        organize_fastq_files_by_directory(positional_arg_count, argc, argv, optind, barcodeFastqFilesString, forwardFastqFilesString, reverseFastqFilesString, &fastq_files, barcode_pattern, forward_pattern, reverse_pattern);
    }
    else{
        organize_fastq_files_by_type(positional_arg_count, argc, argv,optind, barcodeFastqFilesString, forwardFastqFilesString, reverseFastqFilesString, &fastq_files, barcode_pattern, forward_pattern, reverse_pattern,sample_flag);
    }
    whitelist_hash = kh_init(u32ptr);
    read_whiteList(whitelist_filename, whitelist_hash, reverse_complement_whitelist);
    if (sample_barcodes_filename) {
        sample_barcodes = read_features_file(sample_barcodes_filename);
        if (!sample_barcodes) { fprintf(stderr, "Failed to load sample barcodes file %s\n", sample_barcodes_filename); exit(1);}    
        if ((sample_constant_offset_cli == -2) && sample_offset_relative_cli == 0) {
            fprintf(stderr, "Error: must specify --sample_offset or --sample_offset_rel when --sample_barcodes given\n");
            exit(1);
        }
        if (sample_constant_offset_cli != -2 && sample_offset_relative_cli != 0) {
            fprintf(stderr, "Error: specify only ONE of --sample_offset or --sample_offset_rel\n");
            exit(1);
        }
    }
    int demux_nsamples = (sample_barcodes) ? sample_barcodes->number_of_features : 1;
    initialize_unit_sizes();
    const int nSamples=fastq_files.nsamples;
    if (set_search_threads_per_consumer){
        search_threads_per_consumer=set_search_threads_per_consumer;
    }
    if (set_consumer_threads_per_set){
        consumer_threads_per_set=set_consumer_threads_per_set;
    } 
    fprintf(stderr, "Using %d consumer threads and %d search threads per consumer. Max concurrent processes %d\n", consumer_threads_per_set, search_threads_per_consumer, max_concurrent_processes);
    int concurrent_processes=0;
    atomic_int *thread_counter=mmap(NULL, sizeof(atomic_int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    if (thread_counter == MAP_FAILED){
        perror("Failed to allocate memory for thread counter");
        exit(EXIT_FAILURE);
    }
    atomic_init(thread_counter, 0);
    int threads_per_set = consumer_threads_per_set * (1 + search_threads_per_consumer);
    fprintf(stderr,"threads per set %d\n", threads_per_set);
    for (int index=0; index<nSamples; index++){
        const int i=fastq_files.sorted_index[index];
    
        if (concurrent_processes >= max_concurrent_processes){
            wait(NULL);
            concurrent_processes--;
        }
        pid_t pid=fork();
        if (pid< 0){
            perror("Failed to fork");
            exit(EXIT_FAILURE); 
        }
        else if (pid== 0){
            atomic_fetch_add(thread_counter, threads_per_set);
            
            char sample_directory[FILENAME_LENGTH];
            if (sample_flag){
                strcpy(sample_directory, directory);
                strcat(sample_directory, fastq_files.sample_names[i]);
                strcat(sample_directory, "/");
            }
            else{
                strcpy(sample_directory, directory);
                strcat(sample_directory, fastq_files.sample_names[0]);
                strcat(sample_directory, "/");
            }
            fprintf(stderr, "Processing sample directory %s\n", sample_directory);
            if (existing_output_skip(keep_existing, sample_directory)) exit(0);
            sample_args args;
            // The following are now initialized in process_files_in_sample
            // memory_pool_collection *pools=initialize_memory_pool_collection();
            // statistics stats;
            // data_structures hashes;
            // initialize_data_structures(&hashes);
            // initialize_statistics(&stats);
            args.sample_index = i;
            args.directory = sample_directory;
            if (filtered_barcodes_filename)
              args.filtered_barcodes_name = filtered_barcodes_filename;
            else
              args.filtered_barcodes_name = NULL;
            args.fastq_files = &fastq_files;
            args.features = features;
            args.maxHammingDistance = maxHammingDistance;
            args.nThreads = search_threads_per_consumer;
            args.pools = NULL;
            args.stats = NULL;
            args.hashes = NULL;
            args.stringency = stringency;
            args.min_counts = min_counts;
            args.barcode_constant_offset = barcode_constant_offset;
            args.feature_constant_offset = feature_constant_offset;
            args.read_buffer_lines = read_buffer_lines;
            args.average_read_length = average_read_length;
            args.min_posterior = min_posterior;
            args.consumer_threads_per_set = consumer_threads_per_set;
            args.filtered_barcodes_hash = filtered_barcodes_hash;
            args.min_prediction = min_prediction;
            args.min_heatmap = min_heatmap;
            args.demux_nsamples = demux_nsamples;
            args.sample_barcodes = sample_barcodes;
            args.sample_max_hamming = sample_max_hamming;
            args.sample_max_N = sample_max_N;
            args.sample_constant_offset = (sample_constant_offset_cli!=-2)? sample_constant_offset_cli : -1;
            args.sample_offset_relative = sample_offset_relative_cli;
            args.sample_hashes = NULL;
            args.sample_stats = NULL;
            args.sample_pools = NULL;
            process_files_in_sample(&args);
            // cleanup_sample is handled within process_files_in_sample
            atomic_fetch_add(thread_counter, -threads_per_set);
            exit(0);
        }
        concurrent_processes++;

    }
    while (concurrent_processes > 0) {
        wait(NULL);
        concurrent_processes--;
    }
    kh_destroy(u32ptr, whitelist_hash);
    free(whitelist);
    if (barcodeFastqFilesString) free(barcodeFastqFilesString);
    if (forwardFastqFilesString) free(forwardFastqFilesString);
    if (reverseFastqFilesString) free(reverseFastqFilesString);
    if (filtered_barcodes_filename) free(filtered_barcodes_filename);
    if (filtered_barcodes_hash) free_strptr_hash(filtered_barcodes_hash);
    free_feature_arrays(features);
    free_fastq_files_collection(&fastq_files);
    return 0;
}