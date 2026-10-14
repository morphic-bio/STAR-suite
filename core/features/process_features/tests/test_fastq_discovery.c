/* Exercise private triplet discovery without adding a public API. The test
 * linker discards the unused processing routines. */
#include "../src/pf_api.c"
#include <assert.h>
#include <sys/stat.h>
#include <unistd.h>
#include <htslib/bgzf.h>

static const char record[] = "@read\nACGT\n+\nIIII\n";
static void write_read(const char *path, int compression) {
    if (compression == 2) {
        BGZF *f = bgzf_open(path, "w");
        assert(f && bgzf_write(f, record, sizeof(record)-1) == sizeof(record)-1);
        assert(bgzf_close(f) == 0);
    } else if (compression == 1) {
        gzFile f = gzopen(path, "wb");
        assert(f && gzwrite(f, record, sizeof(record)-1) == sizeof(record)-1);
        assert(gzclose(f) == Z_OK);
    } else {
        FILE *f = fopen(path, "w");
        assert(f && fwrite(record, 1, sizeof(record)-1, f) == sizeof(record)-1);
        assert(fclose(f) == 0);
    }
}
static void check_found(const char *dir, int expected) {
    int n = -1;
    char **files = find_files_with_pattern(dir, "_R1_", &n);
    assert(n == expected);
    for (int i = 0; i < n; ++i) {
        char buf[100] = {0};
        gzFile f = gzopen(files[i], "rb");
        assert(f && gzread(f, buf, sizeof(buf)) == sizeof(record)-1);
        assert(strcmp(buf, record) == 0 && gzclose(f) == Z_OK);
    }
    pf_free_file_list(files, n);
}
static void finder(void) {
    char dir[] = "/tmp/pf-discovery-XXXXXX", path[1024], receipt[1024], alias[1024];
    assert(mkdtemp(dir));
    snprintf(receipt, sizeof(receipt), "%s/sample_R1_001.fastq.gz.receipt.json", dir);
    write_read(receipt, 0);
    check_found(dir, 0);
    snprintf(path, sizeof(path), "%s/sample_R1_001.fastq.gz", dir);
    write_read(path, 1);
    check_found(dir, 1);
    assert(unlink(path) == 0 && unlink(receipt) == 0);
    const char *suffixes[] = {".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fastq.bgz",
        ".fq.bgz", ".fastq.bgzf", ".fq.bgzf", ".FASTQ.GZ"};
    for (unsigned i = 0; i < sizeof(suffixes)/sizeof(suffixes[0]); ++i) {
        snprintf(path, sizeof(path), "%s/sample_R1_001%s", dir, suffixes[i]);
        for (int compression = 0; compression <= 2; ++compression) {
            write_read(path, compression);
            check_found(dir, 1);
        }
        assert(unlink(path) == 0);
    }
    /* A FASTQ-named link to a differently named target is still eligible. */
    write_read(receipt, 0);
    snprintf(path, sizeof(path), "%s/link_R1_001.fq", dir);
    assert(symlink(receipt, path) == 0);
    check_found(dir, 1);
    snprintf(alias, sizeof(alias), "%s/link_R1_001.fq.receipt.json", dir);
    assert(symlink(path, alias) == 0);
    check_found(dir, 1);
    assert(unlink(alias) == 0 && unlink(path) == 0 && unlink(receipt) == 0);
    const char *ignored[] = {".txt", ".fastq.gz.tmp", ".fq.bz2", ".fastq.xz", ".gz"};
    for (unsigned i = 0; i < sizeof(ignored)/sizeof(ignored[0]); ++i) {
        snprintf(path, sizeof(path), "%s/sample_R1_001%s", dir, ignored[i]);
        write_read(path, 0);
        check_found(dir, 0);
        assert(unlink(path) == 0);
    }
    snprintf(path, sizeof(path), "%s/directory_R1_001.fastq", dir);
    assert(mkdir(path, 0700) == 0);
    check_found(dir, 0);
    assert(rmdir(path) == 0 && rmdir(dir) == 0);
    puts("PASS finder: supported suffixes and codecs, receipts, symlinks, unsupported names, directories");
}
static void triplets(int explicit_patterns) {
    char dir[] = "/tmp/pf-triplets-XXXXXX", path[1024];
    assert(mkdtemp(dir));
    pf_split_read_layout layout = {0};
    const char *prefix = explicit_patterns ? "custom" : "R";
    if (explicit_patterns) {
        strcpy(layout.fastq_r1_pattern, "_custom1_");
        strcpy(layout.fastq_r2_pattern, "_custom2_");
        strcpy(layout.fastq_r3_pattern, "_custom3_");
    }
    for (int lane = 2; lane >= 1; --lane)
        for (int mate = 1; mate <= 3; ++mate) {
            snprintf(path, sizeof(path), "%s/sample_L00%d_%s%d_001.fastq.gz", dir, lane, prefix, mate);
            write_read(path, 1);
        }
    char **before[3], **after[3];
    int nb = 0, na = 0;
    assert(pf_discover_split_read_triplets(dir, &layout, &before[0], &before[1], &before[2], &nb));
    assert(nb == 2);
    for (int lane = 1; lane <= 2; ++lane)
        for (int mate = 1; mate <= 3; ++mate) {
            snprintf(path, sizeof(path), "%s/sample_L00%d_%s%d_001.fastq.gz.receipt.json", dir, lane, prefix, mate);
            write_read(path, 0);
        }
    assert(pf_discover_split_read_triplets(dir, &layout, &after[0], &after[1], &after[2], &na));
    assert(na == nb);
    if (!explicit_patterns) {
        char **paired[3];
        int counts[3];
        assert(find_paired_fastq_files(dir, "_R1_", "_R2_", "_R3_",
                   &paired[0], &counts[0], &paired[1], &counts[1], &paired[2], &counts[2]) == 0);
        for (int mate = 0; mate < 3; ++mate) {
            assert(counts[mate] == nb);
            for (int lane = 0; lane < nb; ++lane)
                assert(strcmp(paired[mate][lane], before[mate][lane]) == 0);
            pf_free_file_list(paired[mate], counts[mate]);
        }
        puts("PASS barcode-pattern caller: receipt companions excluded");
    }
    for (int mate = 0; mate < 3; ++mate) {
        assert(strcmp(before[mate][0], before[mate][1]) < 0);
        for (int lane = 0; lane < nb; ++lane) {
            assert(strcmp(before[mate][lane], after[mate][lane]) == 0);
            assert(unlink(before[mate][lane]) == 0);
            snprintf(path, sizeof(path), "%s.receipt.json", before[mate][lane]);
            assert(unlink(path) == 0);
        }
        pf_free_file_list(before[mate], nb);
        pf_free_file_list(after[mate], na);
    }
    assert(rmdir(dir) == 0);
    printf("PASS triplets: %s patterns, sorted lanes unchanged with receipts\n", explicit_patterns ? "explicit" : "default");
}
int main(void) { finder(); triplets(0); triplets(1); return 0; }
