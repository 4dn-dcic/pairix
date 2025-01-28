/* The MIT License

   Copyright (c) 2009 Genome Research Ltd (GRL), 2010 Broad Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Soo Lee <duplexa@gmail.com> */

#ifndef __TABIDX_H
#define __TABIDX_H

#define PACKAGE_VERSION "0.3.9"

#include <stdint.h>
#include "kstring.h"
#include "bgzf.h"

#define DEFAULT_REGION_SPLIT_CHARACTER   '|'
#define DEFAULT_REGION_SPLIT_CHARACTER_STR  "|"

#define TI_PRESET_GENERIC 0
#define TI_PRESET_SAM     1
#define TI_PRESET_VCF     2
#define TI_PRESET_PAIRS     3
#define TI_PRESET_MERGED_NODUPS     4
#define TI_PRESET_OLD_MERGED_NODUPS     5

#define TI_FLAG_UCSC      0x10000

typedef int (*ti_fetch_f)(int l, const char *s, void *data);

struct __ti_index_t;
typedef struct __ti_index_t ti_index_t;

struct __ti_iter_t;
typedef struct __ti_iter_t *ti_iter_t;

typedef struct {
	BGZF *fp;
	ti_index_t *idx;
	char *fn, *fnidx;
} pairix_t;

typedef struct {
	int32_t preset;
	int32_t sc, bc, ec; // seq col., beg col. and end col.
	int32_t sc2, bc2, ec2; // seq col., beg col. and end col. for the second coordinate
        char delimiter;  // 2 bytes for 2 chars (this and region_split_character), but 2 bytes padded so that integers are aligned by 4 bytes.
        char region_split_character;
	int32_t meta_char, line_skip;
} ti_conf_t;  // size 40.

typedef struct {
	int beg, end;
	int beg2, end2;
	char *ss, *se;
	char *ss2, *se2;
} ti_interval_t;

typedef struct {
        int tid, beg, end, bin, beg2, end2, bin2;
} ti_intv_t;

typedef struct {
    pairix_t *t;
    ti_iter_t iter;
    int *len;
    char *s;  // This points to iter->str.s. It's redundant but it allows us to flush the string without touching iter->str.s itself.
} iter_unit;

typedef struct {
    iter_unit **iu;
    int n;
    char first;
} merged_iter_t;

typedef struct {
    pairix_t *t;
    ti_iter_t *iter;
    int n;
    int curr;
} sequential_iter_t;


extern ti_conf_t ti_conf_null, ti_conf_gff, ti_conf_bed, ti_conf_psltbl, ti_conf_vcf, ti_conf_sam, ti_conf_pairs, ti_conf_merged_nodups, ti_conf_old_merged_nodups; // preset
extern char global_region_split_character; // separator for 2D (e.g. '|' in "X:1-2000|Y:1-2000")


#ifdef __cplusplus
extern "C" {
#endif

	/*******************
	 * High-level APIs *
	 *******************/

        pairix_t *load_from_file(char *fn);
	pairix_t *ti_open(const char *fn, const char *fnidx);
	int ti_lazy_index_load(pairix_t *t);
	void ti_close(pairix_t *t);
	ti_iter_t ti_query(pairix_t *t, const char *name, int beg, int end);
	sequential_iter_t *ti_query_general(pairix_t *t, const char *name, int beg, int end);
	ti_iter_t ti_queryi(pairix_t *t, int tid, int beg, int end);
	sequential_iter_t *ti_queryi_general(pairix_t *t, int tid, int beg, int end);
	ti_iter_t ti_querys(pairix_t *t, const char *reg);
	sequential_iter_t *ti_querys_general(pairix_t *t, const char *reg);
	ti_iter_t ti_query_2d(pairix_t *t, const char *name, int beg, int end, const char *name2, int beg2, int end2);
	sequential_iter_t *ti_query_2d_general(pairix_t *t, const char *name, int beg, int end, const char *name2, int beg2, int end2);
	ti_iter_t ti_queryi_2d(pairix_t *t, int tid, int beg, int end, int beg2, int end2);
	sequential_iter_t *ti_queryi_2d_general(pairix_t *t, int tid, int beg, int end, int beg2, int end2);
	ti_iter_t ti_querys_2d(pairix_t *t, const char *reg);
        sequential_iter_t *ti_querys_2d_multi(pairix_t *t, const char **regs, int nRegs);
        sequential_iter_t *ti_querys_2d_general(pairix_t *t, const char *reg);

	int ti_query_tid(pairix_t *t, const char *name, int beg, int end);
	int ti_querys_tid(pairix_t *t, const char *reg);
	int ti_query_2d_tid(pairix_t *t, const char *name, int beg, int end, const char *name2, int beg2, int end2);
	int ti_querys_2d_tid(pairix_t *t, const char *reg);
	const char *ti_read(pairix_t *t, ti_iter_t iter, int *len);
        const char *merged_ti_read(merged_iter_t *miter, int *len);
        const char *sequential_ti_read(sequential_iter_t *siter, int *len);

	/* Destroy the iterator */
	void ti_iter_destroy(ti_iter_t iter);
        void destroy_merged_iter(merged_iter_t *miter);

	/* Get the list of sequence names. Each "char*" pointer points to a
	 * internal member of the index, so DO NOT modify the returned
	 * pointer; otherwise the index will be corrupted. The returned
	 * pointer should be freed by a single free() call by the routine
	 * calling this function. The number of sequences is returned at *n. */
	const char **ti_seqname(const ti_index_t *idx, int *n);

        /* get linecount */
        uint64_t get_linecount(const ti_index_t *idx);

        /* get file offset
         * returns number of bgzf blocks spanning a sequence (pair) */
        int get_nblocks(ti_index_t *idx, int tid, BGZF *fp);


        /* check if a pairix-indexed file is a triangle
           ( chromosome pairs occur only in one direction. e.g. if chr1|chr2 exists, chr2|chr1 shouldn't. )
         * returns 1 if triangle
         * returns 0 if not a triangle
         * returns -1 if no chrom (pairs) is found in file
         * returns -2 if the file is 1D-indexed (not applicable) */
        int check_triangle(ti_index_t *idx);


	/******************
	 * Low-level APIs *
	 ******************/

	/* Build the index for file <fn>. File <fn>.px2 will be generated
	 * and overwrite the file of the same name. Return -1 on failure. */
	int ti_index_build(const char *fn, const ti_conf_t *conf);

	/* Load the index from file <fn>.px2. If <fn> is a URL and the index
	 * file is not in the working directory, <fn>.px2 will be
	 * downloaded. Return NULL on failure. */
	ti_index_t *ti_index_load(const char *fn);

	ti_index_t *ti_index_load_local(const char *fnidx);

	/* Destroy the index */
	void ti_index_destroy(ti_index_t *idx);

	/* Parse a region like: chr2, chr2:100, chr2:100-200. Return -1 on failure. */
	int ti_parse_region(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end);
	int ti_parse_region2d(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end, int *begin2, int *end2);

	int ti_get_tid(const ti_index_t *idx, const char *name);

	/* Get the iterator pointing to the first record at the current file
	 * position. If the file is just openned, the iterator points to the
	 * first record in the file. */
	ti_iter_t ti_iter_first(void);

	/* Get the iterator pointing to the first record in region tid:beg-end */
	ti_iter_t ti_iter_query(const ti_index_t *idx, int tid, int beg, int end, int beg2, int end2);

	/* Get the data line pointed by the iterator and iterate to the next record. */
	const char *ti_iter_read(BGZF *fp, ti_iter_t iter, int *len, char seqonly);

	const ti_conf_t *ti_get_conf(ti_index_t *idx);

        /* get column index, 0-based */
        int ti_get_sc(ti_index_t *idx);
        int ti_get_sc2(ti_index_t *idx);
        int ti_get_bc(ti_index_t *idx);
        int ti_get_bc2(ti_index_t *idx);
        int ti_get_ec(ti_index_t *idx);
        int ti_get_ec2(ti_index_t *idx);

        /* get delimiter */
        char ti_get_delimiter(ti_index_t *idx);
        char get_region_split_character(pairix_t *t);
        char ti_get_region_split_character(ti_index_t *idx);

	int ti_get_intv(const ti_conf_t *conf, int len, char *line, ti_interval_t *intv);

        /* convert string 'region1|region2' to 'region2|region1' */
        char* flip_region(char *s, char region_split_character);

        /* create an empty merge_iter_t struct */
        merged_iter_t *create_merged_iter(int n);

        /* fill in an existing (allocated) iter_unit struct from an iter struct */
        void create_iter_unit(pairix_t *t, ti_iter_t iter, iter_unit *iu);

        /* compare two iter_unit structs, compatible with qsort */
        int compare_iter_unit (const void *a, const void *b);

        /* strcmp, argument type modified to be compatible with qsort */
        int strcmp1d(const void* a, const void* b);

        /* double strcmp on the two parts (for string 'xx|yy' vs 'zz|ww' compare 'xx' vs 'zz' first and then 'yy' vs 'ww') */
        int strcmp2d(const void* a, const void* b);

        /* return a uniqified array given an array of strings (generic), returned array must be fried at both array level and element level */
        char **uniq(char** seq_list, int n_seq_list, int *pn_uniq_seq);

        /* given an array of pairix_t structs, get an array of unique key names,
           the returned array must be freed at both array level and element level */
        char** get_unique_merged_seqname(pairix_t **tbs, int n, int *pn_uniq_seq);

        /* get mate1 chromosome list given a list of chromosome pairs, returned array must be freed at both array level and element level. */
        char **get_seq1_list_from_seqpair_list(char** seqpair_list, int n_seqpair_list, int *pn_seq1);

        /* get a sub-list of seq (chrpair) names given seq1, returned array is an array of pointers to the element of the original array */
        char **get_sub_seq_list_for_given_seq1(char *seq1, char **seqpair_list, int n_seqpair_list, int *pn_sub_list);

        /* get a sub-list of seq (chrpair) names given seq2, returned array is an array of pointers to the element of the original array */
        char **get_sub_seq_list_for_given_seq2(char *seq2, char **seqpair_list, int n_seqpair_list, int *pn_sub_list);

        /* get a sub-list of seq2 (chr2) names given seq1, returned array must be freed at both array level and element level. */
        char **get_seq2_list_for_given_seq1(char *seq1, char **seqpair_list, int n_seqpair_list, int *pn_sub_list);

        /* get a sub-list of seq1 (chr1) names given seq2, returned array must be freed at both array level and element level. */
        char **get_seq1_list_for_given_seq2(char *seq2, char **seqpair_list, int n_seqpair_list, int *pn_sub_list);

        /* initialize an empty sequential_iter associated with a pairix_t struct */
        sequential_iter_t *create_sequential_iter(pairix_t *t);

        /* destructor for sequential_iter */
        void destroy_sequential_iter(sequential_iter_t *siter);

        /* add an iter to sequential_iter - the array size is dynamically incremented */
        void add_to_sequential_iter(sequential_iter_t *siter, ti_iter_t iter);

        /* wrapper for 2D string query that allows autoflip */
        sequential_iter_t *querys_2D_wrapper(pairix_t *tb, const char *reg, int flip);


	/*******************
	 * Deprecated APIs *
	 *******************/

	/* The callback version for random access */
	int ti_fetch(BGZF *fp, const ti_index_t *idx, int tid, int beg, int end, void *data, ti_fetch_f func);
	int ti_fetch_2d(BGZF *fp, const ti_index_t *idx, int tid, int beg, int end, int beg2, int end2, void *data, ti_fetch_f func);

	/* Read one line. */
	int ti_readline(BGZF *fp, kstring_t *str);

#ifdef __cplusplus
}
#endif

#endif
