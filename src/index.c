#include <ctype.h>
#include <assert.h>
#include <sys/stat.h>
#include "khash.h"
#include "ksort.h"
#include "kstring.h"
#include "bam_endian.h"
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif
#include "pairix.h"
#include <time.h>

#define TAD_MIN_CHUNK_GAP 32768
#define TAD_LIDX_SHIFT_LARGE_CHR    15
#define TAD_LIDX_SHIFT_ORIGINAL    14
#define TAD_LIDX_SHIFT_PL3 (TAD_LIDX_SHIFT + 3)
#define TAD_LIDX_SHIFT_PL6 (TAD_LIDX_SHIFT + 6)
#define TAD_LIDX_SHIFT_PL9 (TAD_LIDX_SHIFT + 9)
#define TAD_LIDX_SHIFT_PL12 (TAD_LIDX_SHIFT + 12)
#define MAX_CHR_LARGE_CHR 30
#define MAX_CHR_ORIGINAL 29
#define DEFAULT_DELIMITER '\t'
#define MAX_REGION_STR_LEN 10000

int TAD_LIDX_SHIFT = TAD_LIDX_SHIFT_LARGE_CHR;
int MAX_CHR = MAX_CHR_LARGE_CHR;

#define MAGIC_NUMBER "PX2.004\1"
#define OLD_MAGIC_NUMBER2 "PX2.003\1"  // magic number for older version of pairix (0.3.4 - 0.3.5)
#define OLD_MAGIC_NUMBER "PX2.002\1"  // magic number for older version of pairix (up to 0.3.3)


typedef struct {
	uint64_t u, v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)
KSORT_INIT(offt, pair64_t, pair64_lt)

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} ti_binlist_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} ti_lidx_t;

KHASH_MAP_INIT_INT(i, ti_binlist_t)
KHASH_MAP_INIT_STR(s, int)

struct __ti_index_t {
        ti_conf_t conf;
        int32_t n, max;
        khash_t(s) *tname;
        khash_t(i) **index;
        ti_lidx_t *index2;
        uint64_t linecount;
};

struct __ti_iter_t {
        int from_first; // read from the first record; no random access
        int tid, beg, end, beg2, end2, n_off, i, finished;
        uint64_t curr_off;
        kstring_t str;
        const ti_index_t *idx;
        pair64_t *off;
        ti_intv_t intv;
};

ti_conf_t ti_conf_null = { 0, 0, 0, 0, 0, 0, 0, '\t', DEFAULT_REGION_SPLIT_CHARACTER, '#', 0 };
ti_conf_t ti_conf_gff = { 0, 1, 4, 5, 0, 0, 0, '\t', DEFAULT_REGION_SPLIT_CHARACTER, '#', 0 };
ti_conf_t ti_conf_bed = { TI_FLAG_UCSC, 1, 2,  3, 0, 0, 0, '\t', DEFAULT_REGION_SPLIT_CHARACTER, '#', 0 };
ti_conf_t ti_conf_psltbl = { TI_FLAG_UCSC, 15, 17, 18, 0, 0, 0, '\t', DEFAULT_REGION_SPLIT_CHARACTER, '#', 0 };
ti_conf_t ti_conf_sam = { TI_PRESET_SAM, 3, 4, 0, 0, 0, 0, '\t', DEFAULT_REGION_SPLIT_CHARACTER, '@', 0 };
ti_conf_t ti_conf_vcf = { TI_PRESET_VCF, 1, 2, 0, 0, 0, 0, '\t', DEFAULT_REGION_SPLIT_CHARACTER, '#', 0 };
ti_conf_t ti_conf_pairs = { TI_PRESET_PAIRS, 2, 3, 3, 4, 5, 5, '\t', DEFAULT_REGION_SPLIT_CHARACTER, '#', 0 };
ti_conf_t ti_conf_merged_nodups = { TI_PRESET_MERGED_NODUPS, 2, 3, 3, 6, 7, 7, ' ', DEFAULT_REGION_SPLIT_CHARACTER, '#', 0 };
ti_conf_t ti_conf_old_merged_nodups = { TI_PRESET_OLD_MERGED_NODUPS, 3, 4, 4, 7, 8, 8, ' ', DEFAULT_REGION_SPLIT_CHARACTER, '#', 0 };

char global_region_split_character = DEFAULT_REGION_SPLIT_CHARACTER;


/***************
 * read a line *
 ***************/

/*
int ti_readline(BGZF *fp, kstring_t *str)
{
	int c, l = 0;
	str->l = 0;
	while ((c = bgzf_getc(fp)) >= 0 && c != '\n') {
		++l;
		if (c != '\r') kputc(c, str);
	}
	if (c < 0 && l == 0) return -1; // end of file
	return str->l;
}
*/

/* Below is a faster implementation largely equivalent to the one
 * commented out above. */
int ti_readline(BGZF *fp, kstring_t *str)
{
	return bgzf_getline(fp, '\n', str);
}

/*************************************
 * get the interval from a data line *
 *************************************/

static inline int ti_reg2bin(uint32_t beg, uint32_t end)
{
	--end;
	if (beg>>TAD_LIDX_SHIFT == end>>TAD_LIDX_SHIFT) return  4681 + (beg>>TAD_LIDX_SHIFT);
	if (beg>>TAD_LIDX_SHIFT_PL3 == end>>TAD_LIDX_SHIFT_PL3) return   585 + (beg>>TAD_LIDX_SHIFT_PL3);
	if (beg>>TAD_LIDX_SHIFT_PL6 == end>>TAD_LIDX_SHIFT_PL6) return    73 + (beg>>TAD_LIDX_SHIFT_PL6);
	if (beg>>TAD_LIDX_SHIFT_PL9 == end>>TAD_LIDX_SHIFT_PL9) return     9 + (beg>>TAD_LIDX_SHIFT_PL9);
	if (beg>>TAD_LIDX_SHIFT_PL12 == end>>TAD_LIDX_SHIFT_PL12) return     1 + (beg>>TAD_LIDX_SHIFT_PL12);
	return 0;
}

static int get_tid(ti_index_t *idx, const char *ss)
{
	khint_t k;
	int tid;
	k = kh_get(s, idx->tname, ss);
	if (k == kh_end(idx->tname)) { // a new target sequence
		int ret, size;
		// update idx->n, ->max, ->index and ->index2
		if (idx->n == idx->max) {
			idx->max = idx->max? idx->max<<1 : 8;
			idx->index = realloc(idx->index, idx->max * sizeof(void*));
			idx->index2 = realloc(idx->index2, idx->max * sizeof(ti_lidx_t));
		}
		memset(&idx->index2[idx->n], 0, sizeof(ti_lidx_t));
		idx->index[idx->n++] = kh_init(i);
		// update ->tname
		tid = size = kh_size(idx->tname);
		k = kh_put(s, idx->tname, strdup(ss), &ret);
		kh_value(idx->tname, k) = size;
		assert(idx->n == kh_size(idx->tname));
	} else tid = kh_value(idx->tname, k);
	return tid;
}

int ti_get_intv(const ti_conf_t *conf, int len, char *line, ti_interval_t *intv)
{
	int i, b = 0, id = 1, ncols = 0;
	char *s;
	intv->ss = intv->se = 0; intv->ss2 = intv->se2 = 0; intv->beg = intv->end = -1; intv->beg2 = intv->end2 = -1;
	for (i = 0; i <= len; ++i) {
		if (line[i] == conf->delimiter || line[i] == 0) {
                        ++ncols;
			if (id == conf->sc) {
				intv->ss = line + b; intv->se = line + i;
                        } else if (conf->sc2 && id == conf->sc2) {
				intv->ss2 = line + b; intv->se2 = line + i;
			} else if (id == conf->bc) {
				// here ->beg is 0-based.
				intv->beg = intv->end = strtol(line + b, &s, 0);
				if (!(conf->preset&TI_FLAG_UCSC)) --intv->beg;
				else ++intv->end;
				if (intv->beg < 0) intv->beg = 0;
				if (intv->end < 1) intv->end = 1;
			} else if (conf->bc2 && id == conf->bc2) {
				// here ->beg is 0-based.
				intv->beg2 = intv->end2 = strtol(line + b, &s, 0);
				if (!(conf->preset&TI_FLAG_UCSC)) --intv->beg2;
				else ++intv->end2;
				if (intv->beg2 < 0) intv->beg2 = 0;
				if (intv->end2 < 1) intv->end2 = 1;
			} else if(id == conf->ec) {
				if ((conf->preset&0xffff) != TI_PRESET_VCF && (conf->preset&0xffff) != TI_PRESET_SAM) {
					intv->end = strtol(line + b, &s, 0);
				} else if ((conf->preset&0xffff) == TI_PRESET_SAM) {
					if (id == 6) { // CIGAR
						int l = 0, op;
						char *t;
						for (s = line + b; s < line + i;) {
							long x = strtol(s, &t, 10);
							op = toupper(*t);
							if (op == 'M' || op == 'D' || op == 'N') l += x;
							s = t + 1;
						}
						if (l == 0) l = 1;
						intv->end = intv->beg + l;
					}
				} else if ((conf->preset&0xffff) == TI_PRESET_VCF) {
					// FIXME: the following is NOT tested and is likely to be buggy
					if (id == 4) {
						if (b < i) intv->end = intv->beg + (i - b);
					} else if (id == 8) { // look for "END="
						int c = line[i];
						line[i] = 0;
						s = strstr(line + b, "END=");
						if (s == line + b) s += 4;
						else if (s) {
							s = strstr(line + b, ";END=");
							if (s) s += 5;
						}
						if (s) intv->end = strtol(s, &s, 0);
						line[i] = c;
					}
				}
			} else if(conf->ec2 && id == conf->ec2) {
				if ((conf->preset&0xffff) != TI_PRESET_VCF && (conf->preset&0xffff) != TI_PRESET_SAM) {
					intv->end2 = strtol(line + b, &s, 0);
                                }
	      	        } else {
				if ((conf->preset&0xffff) == TI_PRESET_SAM) {
					if (id == 6) { // CIGAR
						int l = 0, op;
						char *t;
						for (s = line + b; s < line + i;) {
							long x = strtol(s, &t, 10);
							op = toupper(*t);
							if (op == 'M' || op == 'D' || op == 'N') l += x;
							s = t + 1;
						}
						if (l == 0) l = 1;
						intv->end = intv->beg + l;
					}
				} else if ((conf->preset&0xffff) == TI_PRESET_VCF) {
					// FIXME: the following is NOT tested and is likely to be buggy
					if (id == 4) {
						if (b < i) intv->end = intv->beg + (i - b);
					} else if (id == 8) { // look for "END="
						int c = line[i];
						line[i] = 0;
						s = strstr(line + b, "END=");
						if (s == line + b) s += 4;
						else if (s) {
							s = strstr(line + b, ";END=");
							if (s) s += 5;
						}
						if (s) intv->end = strtol(s, &s, 0);
						line[i] = c;
					}
				}

                        }
			b = i + 1;
			++id;
		}
	}
	if (intv->ss == 0 || intv->se == 0 || intv->beg < 0 || intv->end < 0 ) return -1;
        if(conf->sc2 && (intv->ss2 == 0 || intv->se2 == 0)) return -1;
        if(conf->bc2 && (intv->beg2 < 0 || intv->end2 < 0)) return -1;
        if(conf->ec2 && (intv->beg2 < 0 || intv->end2 < 0)) return -1;
	return 0;
}


static int get_intv(ti_index_t *idx, kstring_t *str, ti_intv_t *intv)
{
	ti_interval_t x;
        char *str_ptr;
        char sname_double[strlen(str->s)+1];
	intv->tid = intv->beg = intv->end = intv->beg2 = intv->end2 = intv->bin = intv->bin2 =  -1;
        char region_split_character = idx->conf.region_split_character;
	if (ti_get_intv(&idx->conf, str->l, str->s, &x) == 0) {

		char c = *x.se;
                *x.se = '\0';

                if(!x.se2){ //single-chromosome
                  intv->tid = get_tid(idx, x.ss);
                } else { //double-chromosome
		  char c2 = *x.se2;
                  *x.se2 = '\0';
                  strcpy(sname_double,x.ss);
                  str_ptr = sname_double+strlen(sname_double);
                  *str_ptr=region_split_character;
                  str_ptr++;
                  strcpy(str_ptr,x.ss2);
                  intv->tid = get_tid(idx, sname_double);
                  *x.se2=c2;
                }

                *x.se = c;
		intv->beg = x.beg; intv->end = x.end;
		intv->beg2 = x.beg2; intv->end2 = x.end2;
		intv->bin = ti_reg2bin(intv->beg, intv->end);
		intv->bin2 = ti_reg2bin(intv->beg2, intv->end2);

		return (intv->tid >= 0 && intv->beg >= 0 && intv->end >= 0 && ((!idx->conf.bc2 && !idx->conf.ec2) || (intv->beg2 >=0 && intv->end2 >=0))  )? 0 : -1;
	} else {
		fprintf(stderr, "[%s] the following line cannot be parsed and skipped: %s\n", __func__, str->s);
		return -1;
	}
}

/************
 * indexing *
 ************/

// requirement: len <= LEN_MASK
static inline void insert_offset(khash_t(i) *h, int bin, uint64_t beg, uint64_t end)
{
	khint_t k;
	ti_binlist_t *l;
	int ret;
	k = kh_put(i, h, bin, &ret);
	l = &kh_value(h, k);
	if (ret) { // not present
		l->m = 1; l->n = 0;
		l->list = (pair64_t*)calloc(l->m, 16);
	}
	if (l->n == l->m) {
		l->m <<= 1;
		l->list = (pair64_t*)realloc(l->list, l->m * 16);
	}
	l->list[l->n].u = beg; l->list[l->n++].v = end;
}

static inline uint64_t insert_offset2(ti_lidx_t *index2, int _beg, int _end, uint64_t offset)
{
	int i, beg, end;
	beg = _beg >> TAD_LIDX_SHIFT;
	end = (_end - 1) >> TAD_LIDX_SHIFT;
	if (index2->m < end + 1) {
		int old_m = index2->m;
		index2->m = end + 1;
		kroundup32(index2->m);
		index2->offset = (uint64_t*)realloc(index2->offset, index2->m * 8);
		memset(index2->offset + old_m, 0, 8 * (index2->m - old_m));
	}
	if (beg == end) {
		if (index2->offset[beg] == 0) index2->offset[beg] = offset;
	} else {
		for (i = beg; i <= end; ++i)
			if (index2->offset[i] == 0) index2->offset[i] = offset;
	}
	if (index2->n < end + 1) index2->n = end + 1;
	return (uint64_t)beg<<32 | end;
}

static void merge_chunks(ti_index_t *idx)
{
	khash_t(i) *index;
	int i, l, m;
	khint_t k;
	for (i = 0; i < idx->n; ++i) {
		index = idx->index[i];
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			ti_binlist_t *p;
			if (!kh_exist(index, k)) continue;
			p = &kh_value(index, k);
			m = 0;
			for (l = 1; l < p->n; ++l) {
				if (p->list[m].v>>16 == p->list[l].u>>16) p->list[m].v = p->list[l].v;
				else p->list[++m] = p->list[l];
			} // ~for(l)
			p->n = m + 1;
		} // ~for(k)
	} // ~for(i)
}

static void fill_missing(ti_index_t *idx)
{
	int i, j;
	for (i = 0; i < idx->n; ++i) {
		ti_lidx_t *idx2 = &idx->index2[i];
		for (j = 1; j < idx2->n; ++j)
			if (idx2->offset[j] == 0)
				idx2->offset[j] = idx2->offset[j-1];
	}
}

ti_index_t *ti_index_core(BGZF *fp, const ti_conf_t *conf)
{
	int ret;
	ti_index_t *idx;
	uint32_t last_bin, save_bin;
	int32_t last_coor, last_tid, save_tid;
	uint64_t save_off, last_off, lineno = 0, offset0 = (uint64_t)-1, tmp;
	kstring_t *str;

	str = calloc(1, sizeof(kstring_t));

	idx = (ti_index_t*)calloc(1, sizeof(ti_index_t));
	idx->conf = *conf;
	idx->n = idx->max = 0;
	idx->tname = kh_init(s);
	idx->index = 0;
	idx->index2 = 0;
        idx->linecount=0;

	save_bin = save_tid = last_tid = last_bin = 0xffffffffu;
	save_off = last_off = bgzf_tell(fp); last_coor = 0xffffffffu;
	while ((ret = ti_readline(fp, str)) >= 0) {
                idx->linecount++;
		ti_intv_t intv;
		++lineno;
		if (lineno <= idx->conf.line_skip || str->s[0] == idx->conf.meta_char) {
			last_off = bgzf_tell(fp);
			continue;
		}
		get_intv(idx, str, &intv);
                if ( intv.beg<0 || intv.end<0 )
                {
                    fprintf(stderr,"[ti_index_core] the indexes overlap or are out of bounds\n");
                    return(NULL);
                }
		if (last_tid != intv.tid) { // change of chromosomes
                if (last_tid>intv.tid )
                {
                    fprintf(stderr,"[ti_index_core] the chromosome blocks not continuous at line %llu, is the file sorted? [pos %d]\n",(unsigned long long)lineno,intv.beg+1);
                    return(NULL);
                }
		    last_tid = intv.tid;
		    last_bin = 0xffffffffu;
		} else if (last_coor > intv.beg) {
		    fprintf(stderr, "[ti_index_core] the file out of order at line %llu\n", (unsigned long long)lineno);
		    return(NULL);
		}
		tmp = insert_offset2(&idx->index2[intv.tid], intv.beg, intv.end, last_off);
		if (last_off == 0) offset0 = tmp;
		if (intv.bin != last_bin) { // then possibly write the binning index
			if (save_bin != 0xffffffffu) // save_bin==0xffffffffu only happens to the first record
				insert_offset(idx->index[save_tid], save_bin, save_off, last_off);
			save_off = last_off;
			save_bin = last_bin = intv.bin;
			save_tid = intv.tid;
			if (save_tid < 0) break;
		}
		if (bgzf_tell(fp) <= last_off) {
			fprintf(stderr, "[ti_index_core] bug in BGZF: %llx < %llx\n",
					(unsigned long long)bgzf_tell(fp), (unsigned long long)last_off);
			return(NULL);
		}
		last_off = bgzf_tell(fp);
		last_coor = intv.beg;
	}
	if (save_tid >= 0) insert_offset(idx->index[save_tid], save_bin, save_off, bgzf_tell(fp));
	merge_chunks(idx);
	fill_missing(idx);
	if (offset0 != (uint64_t)-1 && idx->n && idx->index2[0].offset) {
		int i, beg = offset0>>32, end = offset0&0xffffffffu;
		for (i = beg; i <= end; ++i) idx->index2[0].offset[i] = 0;
	}

	free(str->s); free(str);
	return idx;
}

void ti_index_destroy(ti_index_t *idx)
{
	khint_t k;
	int i;
	if (idx == 0) return;
	// destroy the name hash table
	for (k = kh_begin(idx->tname); k != kh_end(idx->tname); ++k) {
		if (kh_exist(idx->tname, k))
			free((char*)kh_key(idx->tname, k));
	}
	kh_destroy(s, idx->tname);
	// destroy the binning index
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index = idx->index[i];
		ti_lidx_t *index2 = idx->index2 + i;
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			if (kh_exist(index, k))
				free(kh_value(index, k).list);
		}
		kh_destroy(i, index);
		free(index2->offset);
	}
	free(idx->index);
	// destroy the linear index
	free(idx->index2);
	free(idx);
}

/******************
 * index file I/O *
 ******************/

void ti_index_save(const ti_index_t *idx, BGZF *fp)
{
	int32_t i, size, ti_is_be;
	khint_t k;
	ti_is_be = bam_is_big_endian();
	bgzf_write(fp, MAGIC_NUMBER, 8);
	if (ti_is_be) {
		uint32_t x = idx->n;
		bgzf_write(fp, bam_swap_endian_4p(&x), 4);
	} else bgzf_write(fp, &idx->n, 4);
	if (ti_is_be) {
		uint64_t x = idx->linecount;
		bgzf_write(fp, bam_swap_endian_8p(&x), 8);
	} else bgzf_write(fp, &idx->linecount, 8);
	assert(sizeof(ti_conf_t) == 40);
	if (ti_is_be) { // write ti_conf_t;
		uint32_t x[6];
		memcpy(x, &idx->conf, 40);
		for (i = 0; i < 6; ++i) bgzf_write(fp, bam_swap_endian_4p(&x[i]), 4);
	} else bgzf_write(fp, &idx->conf, sizeof(ti_conf_t));
	{ // write target names
		char **name;
		int32_t l = 0;
		name = calloc(kh_size(idx->tname), sizeof(void*));
		for (k = kh_begin(idx->tname); k != kh_end(idx->tname); ++k)
			if (kh_exist(idx->tname, k))
				name[kh_value(idx->tname, k)] = (char*)kh_key(idx->tname, k);
		for (i = 0; i < kh_size(idx->tname); ++i)
			l += strlen(name[i]) + 1;
		if (ti_is_be) bgzf_write(fp, bam_swap_endian_4p(&l), 4);
		else bgzf_write(fp, &l, 4);
		for (i = 0; i < kh_size(idx->tname); ++i)
			bgzf_write(fp, name[i], strlen(name[i]) + 1);
		free(name);
	}
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index = idx->index[i];
		ti_lidx_t *index2 = idx->index2 + i;
		// write binning index
		size = kh_size(index);
		if (ti_is_be) { // big endian
			uint32_t x = size;
			bgzf_write(fp, bam_swap_endian_4p(&x), 4);
		} else bgzf_write(fp, &size, 4);
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			if (kh_exist(index, k)) {
				ti_binlist_t *p = &kh_value(index, k);
				if (ti_is_be) { // big endian
					uint32_t x;
					x = kh_key(index, k); bgzf_write(fp, bam_swap_endian_4p(&x), 4);
					x = p->n; bgzf_write(fp, bam_swap_endian_4p(&x), 4);
					for (x = 0; (int)x < p->n; ++x) {
						bam_swap_endian_8p(&p->list[x].u);
						bam_swap_endian_8p(&p->list[x].v);
					}
					bgzf_write(fp, p->list, 16 * p->n);
					for (x = 0; (int)x < p->n; ++x) {
						bam_swap_endian_8p(&p->list[x].u);
						bam_swap_endian_8p(&p->list[x].v);
					}
				} else {
					bgzf_write(fp, &kh_key(index, k), 4);
					bgzf_write(fp, &p->n, 4);
					bgzf_write(fp, p->list, 16 * p->n);
				}
			}
		}
		// write linear index (index2)
		if (ti_is_be) {
			int x = index2->n;
			bgzf_write(fp, bam_swap_endian_4p(&x), 4);
		} else bgzf_write(fp, &index2->n, 4);
		if (ti_is_be) { // big endian
			int x;
			for (x = 0; (int)x < index2->n; ++x)
				bam_swap_endian_8p(&index2->offset[x]);
			bgzf_write(fp, index2->offset, 8 * index2->n);
			for (x = 0; (int)x < index2->n; ++x)
				bam_swap_endian_8p(&index2->offset[x]);
		} else bgzf_write(fp, index2->offset, 8 * index2->n);
	}
}

static ti_index_t *ti_index_load_core(BGZF *fp)
{
	int i, ti_is_be;
	char magic[8];
	ti_index_t *idx;
	ti_is_be = bam_is_big_endian();
	if (fp == 0) {
		fprintf(stderr, "[ti_index_load_core] fail to load index.\n");
		return 0;
	}
	bgzf_read(fp, magic, 8);
	if (strncmp(magic, MAGIC_NUMBER, 8)) {
            if (strncmp(magic, OLD_MAGIC_NUMBER, 8)==0) {
                TAD_LIDX_SHIFT = TAD_LIDX_SHIFT_ORIGINAL;
                MAX_CHR = MAX_CHR_ORIGINAL; 
            }
            else if(strncmp(magic, OLD_MAGIC_NUMBER2, 8)==0) {
            }
            else {
		fprintf(stderr, "[ti_index_load] wrong magic number. Re-index if your index file was created by an earlier version of pairix.\n");
		return 0;
            }
	}
	idx = (ti_index_t*)calloc(1, sizeof(ti_index_t));
	bgzf_read(fp, &idx->n, 4);
	if (ti_is_be) bam_swap_endian_4p(&idx->n);
        if(strncmp(magic, MAGIC_NUMBER, 8)==0) {
      	    bgzf_read(fp, &idx->linecount, 8);
	    if (ti_is_be) bam_swap_endian_8p(&idx->linecount);
        }
        else if(strncmp(magic, OLD_MAGIC_NUMBER2, 8)==0 || strncmp(magic, OLD_MAGIC_NUMBER, 8)==0) {
            bgzf_read(fp, &idx->linecount, 4);
            if (ti_is_be) bam_swap_endian_4p(&idx->linecount);
        }
	idx->tname = kh_init(s);
	idx->index = (khash_t(i)**)calloc(idx->n, sizeof(void*));
	idx->index2 = (ti_lidx_t*)calloc(idx->n, sizeof(ti_lidx_t));
	// read idx->conf
	bgzf_read(fp, &idx->conf, sizeof(ti_conf_t));
	if (ti_is_be) {
		bam_swap_endian_4p(&idx->conf.preset);
		bam_swap_endian_4p(&idx->conf.sc);
		bam_swap_endian_4p(&idx->conf.bc);
		bam_swap_endian_4p(&idx->conf.ec);
		bam_swap_endian_4p(&idx->conf.sc2);
		bam_swap_endian_4p(&idx->conf.bc2);
		bam_swap_endian_4p(&idx->conf.ec2);
		bam_swap_endian_4p(&idx->conf.delimiter);
		bam_swap_endian_4p(&idx->conf.meta_char);
		bam_swap_endian_4p(&idx->conf.line_skip);
	}
	{ // read target names
		int j, ret;
		kstring_t *str;
		int32_t l;
		uint8_t *buf;
		bgzf_read(fp, &l, 4);
		if (ti_is_be) bam_swap_endian_4p(&l);
		buf = calloc(l, 1);
		bgzf_read(fp, buf, l);
		str = calloc(1, sizeof(kstring_t));
		for (i = j = 0; i < l; ++i) {
			if (buf[i] == 0) {
				khint_t k = kh_put(s, idx->tname, strdup(str->s), &ret);
				kh_value(idx->tname, k) = j++;
				str->l = 0;
			} else kputc(buf[i], str);
		}
		free(str->s); free(str); free(buf);
	}
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index;
		ti_lidx_t *index2 = idx->index2 + i;
		uint32_t key, size;
		khint_t k;
		int j, ret;
		ti_binlist_t *p;
		index = idx->index[i] = kh_init(i);
		// load binning index
		bgzf_read(fp, &size, 4);
		if (ti_is_be) bam_swap_endian_4p(&size);
		for (j = 0; j < (int)size; ++j) {
			bgzf_read(fp, &key, 4);
			if (ti_is_be) bam_swap_endian_4p(&key);
			k = kh_put(i, index, key, &ret);
			p = &kh_value(index, k);
			bgzf_read(fp, &p->n, 4);
			if (ti_is_be) bam_swap_endian_4p(&p->n);
			p->m = p->n;
			p->list = (pair64_t*)malloc(p->m * 16);
			bgzf_read(fp, p->list, 16 * p->n);
			if (ti_is_be) {
				int x;
				for (x = 0; x < p->n; ++x) {
					bam_swap_endian_8p(&p->list[x].u);
					bam_swap_endian_8p(&p->list[x].v);
				}
			}
		}
		// load linear index
		bgzf_read(fp, &index2->n, 4);
		if (ti_is_be) bam_swap_endian_4p(&index2->n);
		index2->m = index2->n;
		index2->offset = (uint64_t*)calloc(index2->m, 8);
		bgzf_read(fp, index2->offset, index2->n * 8);
		if (ti_is_be)
			for (j = 0; j < index2->n; ++j) bam_swap_endian_8p(&index2->offset[j]);
	}
	return idx;
}

ti_index_t *ti_index_load_local(const char *fnidx)
{
	BGZF *fp;
	fp = bgzf_open(fnidx, "r");
	if (fp) {
		ti_index_t *idx = ti_index_load_core(fp);
		bgzf_close(fp);
		return idx;
	} else return 0;
}

#ifdef _USE_KNETFILE
static void download_from_remote(const char *url)
{
	const int buf_size = 1 * 1024 * 1024;
	char *fn;
	FILE *fp;
	uint8_t *buf;
	knetFile *fp_remote;
	int l;
	if (strstr(url, "ftp://") != url && strstr(url, "http://") != url) return;
	l = strlen(url);
	for (fn = (char*)url + l - 1; fn >= url; --fn)
		if (*fn == '/') break;
	++fn; // fn now points to the file name
	fp_remote = knet_open(url, "r");
	if (fp_remote == 0) {
		fprintf(stderr, "[download_from_remote] fail to open remote file.\n");
		return;
	}
	if ((fp = fopen(fn, "w")) == 0) {
		fprintf(stderr, "[download_from_remote] fail to create file in the working directory.\n");
		knet_close(fp_remote);
		return;
	}
	buf = (uint8_t*)calloc(buf_size, 1);
	while ((l = knet_read(fp_remote, buf, buf_size)) != 0)
		fwrite(buf, 1, l, fp);
	free(buf);
	fclose(fp);
	knet_close(fp_remote);
}
#else
static void download_from_remote(const char *url)
{
	return;
}
#endif

static char *get_local_version(const char *fn)
{
        struct stat sbuf;
	char *fnidx = (char*)calloc(strlen(fn) + 5, 1);
	strcat(strcpy(fnidx, fn), ".px2");
	if ((strstr(fnidx, "ftp://") == fnidx || strstr(fnidx, "http://") == fnidx)) {
		char *p, *url;
		int l = strlen(fnidx);
		for (p = fnidx + l - 1; p >= fnidx; --p)
			if (*p == '/') break;
		url = fnidx; fnidx = strdup(p + 1);
		if (stat(fnidx, &sbuf) == 0) {
			free(url);
			return fnidx;
		}
		fprintf(stderr, "[%s] downloading the index file...\n", __func__);
		download_from_remote(url);
		free(url);
	}
        if (stat(fnidx, &sbuf) == 0) return fnidx;
	free(fnidx);
        return 0;
}

const char **ti_seqname(const ti_index_t *idx, int *n)
{
	const char **names;
	khint_t k;
	*n = idx->n;
	names = calloc(idx->n, sizeof(void*));
	for (k = kh_begin(idx->tname); k < kh_end(idx->tname); ++k)
		if (kh_exist(idx->tname, k))
			names[kh_val(idx->tname, k)] = kh_key(idx->tname, k);
	return names;
}

ti_index_t *ti_index_load(const char *fn)
{
	ti_index_t *idx;
        char *fname = get_local_version(fn);
	if (fname == 0) return 0;
	idx = ti_index_load_local(fname);
	if (idx == 0) fprintf(stderr, "[ti_index_load] fail to load the index: %s\n", fname);
        if(fname) free(fname);
	return idx;
}

int ti_index_build2(const char *fn, const ti_conf_t *conf, const char *_fnidx)
{
	char *fnidx;
	BGZF *fp, *fpidx;
	ti_index_t *idx;
	if ((fp = bgzf_open(fn, "r")) == 0) {
		fprintf(stderr, "[ti_index_build2] fail to open the file: %s\n", fn);
		return -1;
	}
	idx = ti_index_core(fp, conf);
        if(!idx) return -1;
	bgzf_close(fp);
	if (_fnidx == 0) {
		fnidx = (char*)calloc(strlen(fn) + 5, 1);
		strcpy(fnidx, fn); strcat(fnidx, ".px2");
	} else fnidx = strdup(_fnidx);
	fpidx = bgzf_open(fnidx, "w");
	if (fpidx == 0) {
		fprintf(stderr, "[ti_index_build2] fail to create the index file.\n");
		free(fnidx);
		return -1;
	}
	ti_index_save(idx, fpidx);
	ti_index_destroy(idx);
	bgzf_close(fpidx);
	free(fnidx);
	return 0;
}

int ti_index_build(const char *fn, const ti_conf_t *conf)
{
	return ti_index_build2(fn, conf, 0);
}

/********************************************
 * parse a region in the format chr:beg-end *
 ********************************************/

int ti_get_tid(const ti_index_t *idx, const char *name)
{
	khiter_t iter;
	const khash_t(s) *h = idx->tname;
	iter = kh_get(s, h, name); /* get the tid */
	if (iter == kh_end(h)) return -1;
	return kh_value(h, iter);
}

int ti_parse_region(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end)
{
	char *s, *p;
	int i, l, k;
	l = strlen(str);
	p = s = (char*)malloc(l+1);
	/* squeeze out "," */
	for (i = k = 0; i != l; ++i)
		if (str[i] != ',' && !isspace(str[i])) s[k++] = str[i];
	s[k] = 0;
	for (i = 0; i != k; ++i) if (s[i] == ':') break;
	s[i] = 0;
	if ((*tid = ti_get_tid(idx, s)) < 0) {
		free(s);
		return -1;
	}
	if (i == k) { /* dump the whole sequence */
		*begin = 0; *end = 1<<MAX_CHR; free(s);
		return 0;
	}
	for (p = s + i + 1; i != k; ++i) if (s[i] == '-') break;
	*begin = atoi(p);
	if (i < k) {
		p = s + i + 1;
		*end = atoi(p);
	} else *end = 1<<MAX_CHR;
	if (*begin > 0) --*begin;
	free(s);
	if (*begin > *end) return -2;

	return 0;
}


// thie function can handle both 1d and 2d query automatically
// if 1d, begin2 and end2 will have value -1.
// query string error: -1
// region_split_character not matching error: -2
// memory allocation error: -3
int ti_parse_region2d(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end, int *begin2, int *end2)
{
	char *s, *p, *sname, *tmp_s;
	int i, l, k, h;
        int coord1s, coord1e, coord2s, coord2e, pos1s, pos2s;
        char region_split_character = ti_get_region_split_character(idx);

	l = strlen(str);
	p = s = (char*)malloc(l+1);
	/* squeeze out "," */
	for (i = k = 0; i != l; ++i)
	  if (str[i] != ',' && !isspace(str[i])) s[k++] = str[i];
	s[k] = 0;

        /* split by dimension */
        for(i = 0; i != k; i++) if( s[i] == region_split_character) break;
        s[i]=0;

        /* get data dimension */
        int dim = ti_get_sc2(idx)+1==0?1:2;

        if(i == k && dim == 1) { //1d query
          *begin2=-1; *end2=-1;
          int res = ti_parse_region(idx,str,tid,begin,end);
          free(s); return (res);
        }
        if(i == k && dim == 2) { //1d query on 2d data : interprete query 'x' as 'x|x'
          tmp_s = (char*)realloc(s, k*2+2);
          if(tmp_s) s = tmp_s;
          else return(-3);  // memory alloc error
          strcpy(s+i+1, s);
          s[i] = region_split_character;
          k = k*2+1;
        }

        //2d query on 2d data
        coord1s=0; coord1e=i; coord2s=i+1; coord2e=k;

        /* split into chr and pos */
        for (i = coord1s; i != coord1e; ++i) if (s[i] == ':') break;
        s[i]=0; pos1s=i+1;
        for (i = coord2s; i != coord2e; ++i) if (s[i] == ':') break;
        s[i]=0; pos2s=i+1;

        /* concatenate chromosomes */
        sname = (char*)malloc(l+1);
        strcpy(sname, s + coord1s);
        h=strlen(sname);
        sname[h]= region_split_character;
        strcpy(sname+h+1, s+coord2s);


        if ((*tid = ti_get_tid(idx, sname)) < 0) {
	    free(s); free(sname);
	    return -1;
        }

        /* parsing pos1 */
        if (pos1s-1 == coord1e) { /* dump the whole sequence */
	    *begin = 0; *end = 1<<MAX_CHR;
        } else {
            p = s + pos1s;
	    for (i = pos1s ; i != coord1e; ++i) if (s[i] == '-') break;
	    *begin = atoi(p);
	    if (i < coord1e) {
		p = s + i + 1;
		*end = atoi(p);
  	    } else *end = 1<<MAX_CHR;
  	    if (*begin > 0) --*begin;
        }

        /* parsing pos2 */
        if (pos2s-1 == coord2e) { /* dump the whole sequence */
		*begin2 = 0; *end2 = 1<<MAX_CHR;
        } else {
            p = s + pos2s;
	    for (i = pos2s ; i != coord2e; ++i) if (s[i] == '-') break;
	    *begin2 = atoi(p);
	    if (i < coord2e) {
		p = s + i + 1;
		*end2 = atoi(p);
	    } else *end2 = 1<<MAX_CHR;
	    if (*begin2 > 0) --*begin2;
        }

	free(s); free(sname);
	if (*begin > *end) return -1;
   	if (*begin2!=-1 && *begin2 > *end2) return -1;
 	return 0;
}




/*******************************
 * retrieve a specified region *
 *******************************/

#define MAX_BIN 37450 // =(8^6-1)/7+1

// #define MAX_BIN 74898
// #define MAX_BIN 149794
// #define MAX_BIN 299594

static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[MAX_BIN])
{
	int i = 0, k;
	if (beg >= end) return 0;
	if (end > 1u<<MAX_CHR) {
            end = 1u<<MAX_CHR;
            fprintf(stderr, "Warning: maximum chromosome size is 2^%d.\n", MAX_CHR);
            if(MAX_CHR == MAX_CHR_ORIGINAL) fprintf(stderr, "Old version of index detected. Re-index to increase the chromosomze size limit to 2^%d.\n", MAX_CHR_LARGE_CHR);
        }
	--end;
	list[i++] = 0;
	for (k =     1 + (beg>>TAD_LIDX_SHIFT_PL12); k <=     1 + (end>>TAD_LIDX_SHIFT_PL12); ++k) list[i++] = k;
	for (k =     9 + (beg>>TAD_LIDX_SHIFT_PL9); k <=     9 + (end>>TAD_LIDX_SHIFT_PL9); ++k) list[i++] = k;
	for (k =    73 + (beg>>TAD_LIDX_SHIFT_PL6); k <=    73 + (end>>TAD_LIDX_SHIFT_PL6); ++k) list[i++] = k;
	for (k =   585 + (beg>>TAD_LIDX_SHIFT_PL3); k <=   585 + (end>>TAD_LIDX_SHIFT_PL3); ++k) list[i++] = k;
	for (k =  4681 + (beg>>TAD_LIDX_SHIFT); k <=  4681 + (end>>TAD_LIDX_SHIFT); ++k) list[i++] = k;
	return i;
}

ti_iter_t ti_iter_first()
{
	ti_iter_t iter;
	iter = calloc(1, sizeof(struct __ti_iter_t));
	iter->from_first = 1;
	return iter;
}


uint64_t get_linecount(const ti_index_t *idx)
{
        return(idx->linecount);
}


ti_iter_t ti_iter_query(const ti_index_t *idx, int tid, int beg, int end, int beg2, int end2 ){ //beg2, end2 should be -1 for 1d query.
	uint16_t *bins;
	int i, n_bins, n_off;
	pair64_t *off;
	khint_t k;
	khash_t(i) *index;
	uint64_t min_off;
	ti_iter_t iter = 0;

	if (beg < 0) beg = 0;
	if (end < beg) return 0;
	// initialize the iterator
	iter = calloc(1, sizeof(struct __ti_iter_t));
	iter->idx = idx; iter->tid = tid; iter->beg = beg; iter->end = end; iter->beg2 = beg2; iter->end2 = end2;  iter->i = -1;
	// random access
	bins = (uint16_t*)calloc(MAX_BIN, 2);
	n_bins = reg2bins(beg, end, bins);
	index = idx->index[tid];
	if (idx->index2[tid].n > 0) {
		min_off = (beg>>TAD_LIDX_SHIFT >= idx->index2[tid].n)? idx->index2[tid].offset[idx->index2[tid].n-1]
			: idx->index2[tid].offset[beg>>TAD_LIDX_SHIFT];
		if (min_off == 0) { // improvement for index files built by tabix prior to 0.1.4
			int n = beg>>TAD_LIDX_SHIFT;
			if (n > idx->index2[tid].n) n = idx->index2[tid].n;
			for (i = n - 1; i >= 0; --i)
				if (idx->index2[tid].offset[i] != 0) break;
			if (i >= 0) min_off = idx->index2[tid].offset[i];
		}
	} else min_off = 0; // tabix 0.1.2 may produce such index files
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(i, index, bins[i])) != kh_end(index))
			n_off += kh_value(index, k).n;
	}
	if (n_off == 0) {
		free(bins); return iter;
	}
	off = (pair64_t*)calloc(n_off, 16);
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(i, index, bins[i])) != kh_end(index)) {
			int j;
			ti_binlist_t *p = &kh_value(index, k);
			for (j = 0; j < p->n; ++j)
				if (p->list[j].v > min_off) off[n_off++] = p->list[j];
		}
	}
	if (n_off == 0) {
		free(bins); free(off); return iter;
	}
	free(bins);
	{
		int l;
		ks_introsort(offt, n_off, off);
		// resolve completely contained adjacent blocks
		for (i = 1, l = 0; i < n_off; ++i)
			if (off[l].v < off[i].v)
				off[++l] = off[i];
		n_off = l + 1;
		// resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
		for (i = 1; i < n_off; ++i)
			if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
		{ // merge adjacent blocks
			for (i = 1, l = 0; i < n_off; ++i) {
				if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
				else off[++l] = off[i];
			}
			n_off = l + 1;
		}
	}
	iter->n_off = n_off; iter->off = off;
	return iter;
}


const char *ti_iter_read(BGZF *fp, ti_iter_t iter, int *len, char seqonly)
{
        if (!iter) return 0;
	if (iter->finished) return 0;
	if (iter->from_first) {
		int ret;
		if ((ret = ti_readline(fp, &iter->str)) < 0) {
			iter->finished = 1;
			return 0;
		} else {
			if (len) *len = iter->str.l;
			return iter->str.s;
		}
	}
	if (iter->n_off == 0) return 0;
	while (1) {
		int ret;
		if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
			if (iter->i == iter->n_off - 1) break; // no more chunks
			if (iter->i >= 0) assert(iter->curr_off == iter->off[iter->i].v); // otherwise bug
			if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
				bgzf_seek(fp, iter->off[iter->i+1].u, SEEK_SET);
				iter->curr_off = bgzf_tell(fp);
			}
			++iter->i;
		}
		if ((ret = ti_readline(fp, &iter->str)) >= 0) {
			iter->curr_off = bgzf_tell(fp);
			if (iter->str.s[0] == iter->idx->conf.meta_char) continue;
			get_intv((ti_index_t*)iter->idx, &iter->str, &iter->intv);
                        if(seqonly)
                                if(iter->intv.tid == iter->tid) {
                                      if (len) *len = iter->str.l;
                                      return iter->str.s;  // compare only chromosome (chromosome pair) not position.
                                } else break;
			else if (iter->intv.tid != iter->tid || iter->intv.beg >= iter->end ) break; // no need to proceed
			else if (iter->intv.end > iter->beg && iter->end > iter->intv.beg && ( iter->beg2==-1 || iter->end2==-1 || (iter->intv.end2 > iter->beg2 && iter->end2 > iter->intv.beg2) )) {
				if (len) *len = iter->str.l;
				return iter->str.s;
			} else continue;
		} else break; // end of file
	}
	iter->finished = 1;
	return 0;
}

void ti_iter_destroy(ti_iter_t iter)
{
	if (iter) {
		free(iter->str.s); free(iter->off);
		free(iter);
	}
}

int ti_fetch(BGZF *fp, const ti_index_t *idx, int tid, int beg, int end, void *data, ti_fetch_f func)
{
	ti_iter_t iter;
	const char *s;
	int len;
	iter = ti_iter_query(idx, tid, beg, end, -1, -1);
	while ((s = ti_iter_read(fp, iter, &len, 0)) != 0)
		func(len, s, data);
	ti_iter_destroy(iter);
	return 0;
}

int ti_fetch_2d(BGZF *fp, const ti_index_t *idx, int tid, int beg, int end, int beg2, int end2, void *data, ti_fetch_f func)
{
	ti_iter_t iter;
	const char *s;
	int len;
	iter = ti_iter_query(idx, tid, beg, end, beg2, end2);
	while ((s = ti_iter_read(fp, iter, &len, 0)) != 0)
		func(len, s, data);
	ti_iter_destroy(iter);
	return 0;
}

const ti_conf_t *ti_get_conf(ti_index_t *idx) { return idx? &idx->conf : 0; }

/* Get 0-based column index and delimiter */
int ti_get_sc(ti_index_t *idx) { return idx? idx->conf.sc-1 : -1; }
int ti_get_sc2(ti_index_t *idx) { return idx? idx->conf.sc2-1 : -1; }
int ti_get_bc(ti_index_t *idx) { return idx? idx->conf.bc-1 : -1; }
int ti_get_bc2(ti_index_t *idx) { return idx? idx->conf.bc2-1 : -1; }
int ti_get_ec(ti_index_t *idx) { return idx? idx->conf.ec-1 : -1; }
int ti_get_ec2(ti_index_t *idx) { return idx? idx->conf.ec2-1 : -1; }
char ti_get_delimiter(ti_index_t *idx) { return idx? idx->conf.delimiter : 0; }
char ti_get_region_split_character(ti_index_t *idx){ return idx? idx->conf.region_split_character : 0; }


/*******************
 * High-level APIs *
 *******************/

pairix_t *ti_open(const char *fn, const char *fnidx)
{
	pairix_t *t;
	BGZF *fp;
	if ((fp = bgzf_open(fn, "r")) == 0) return 0;
	t = calloc(1, sizeof(pairix_t));
	t->fn = strdup(fn);
	if (fnidx) t->fnidx = strdup(fnidx);
	t->fp = fp;
	return t;
}

void ti_close(pairix_t *t)
{
	if (t) {
		bgzf_close(t->fp);
		if (t->idx) ti_index_destroy(t->idx);
		free(t->fn); free(t->fnidx);
		free(t);
	}
}

int ti_lazy_index_load(pairix_t *t)
{
	if (t->idx == 0) { // load index
		if (t->fnidx) t->idx = ti_index_load_local(t->fnidx);
		else t->idx = ti_index_load(t->fn);
		if (t->idx == 0) return -1; // fail to load index
	}
	return 0;
}


// reg can contain wildcard for one of the two mates ('*')
sequential_iter_t *ti_querys_2d_general(pairix_t *t, const char *reg)
{
   int n_seqpair_list;
   int n_sub_list;
   char *sp, **chr1list, **chr2list;
   char *chrend;
   char chronly=1;
   int i;

   char region_split_character = t->idx->conf.region_split_character;

   if((sp = strchr(reg, region_split_character)) != NULL){
      if(sp == reg + 1 && reg[0]=='*') {    // '*|c:s-e'
         char *chr2 = sp + 1;

         // get the second chromosome and the list of the first chromosomes pairing with the second chromosome.
         // extract only chr part temporarily (split the beg and end part for now)
         if( (chrend = strchr(chr2, ':')) != NULL) {
            *chrend=0; chronly=0;
         }
         char **chrpairlist = (char**)ti_seqname(t->idx, &n_seqpair_list);
         chr1list = get_seq1_list_for_given_seq2(chr2, chrpairlist, n_seqpair_list, &n_sub_list);
         if(chronly==0) *chrend=':';  // revert to original region string including beg and end
         // create an array of regions in string.
         char **regions = malloc(n_sub_list * sizeof(char*));
         for(i=0;i<n_sub_list;i++){
            regions[i] = malloc((strlen(chr1list[i]) + strlen(chr2) + 2) * sizeof(char));
            strcpy(regions[i], chr1list[i]);
            *(regions[i] + strlen(regions[i]) + 1) = 0;
            *(regions[i] + strlen(regions[i])) = region_split_character;
            strcat(regions[i], chr2);
         }
         free(chrpairlist);
         for(i=0;i<n_sub_list;i++) free(chr1list[i]);
         free(chr1list);

         // multi-region query
         sequential_iter_t *siter = ti_querys_2d_multi(t, (const char**)regions, n_sub_list);
         for(i=0;i<n_sub_list;i++) free(regions[i]);
         free(regions);
         return(siter);

      } else if(strlen(sp) == 2 && sp[1]=='*'){  // 'c:s-e|*'
         *sp=0; char *chr1 = reg;
         if( (chrend = strchr(chr1, ':')) != NULL) {
            *chrend=0; chronly=0;
         }
         char **chrpairlist = (char**)ti_seqname(t->idx, &n_seqpair_list);
         chr2list = get_seq2_list_for_given_seq1(chr1, chrpairlist, n_seqpair_list, &n_sub_list);
         if(chronly==0) *chrend=':';  // revert to original region string including beg and end

         // create an array of regions in string.
         char **regions = malloc(n_sub_list * sizeof(char*));
         for(i=0;i<n_sub_list;i++){
            regions[i] = malloc((strlen(chr2list[i]) + strlen(chr1) + 2) * sizeof(char));
            strcpy(regions[i], chr1);
            *(regions[i] + strlen(regions[i]) + 1) = 0;
            *(regions[i] + strlen(regions[i])) = region_split_character;
            strcat(regions[i], chr2list[i]);
         }
         free(chrpairlist);
         for(i=0;i<n_sub_list;i++) free(chr2list[i]);
         free(chr2list);

         sequential_iter_t *siter = ti_querys_2d_multi(t, (const char**)regions, n_sub_list);
         for(i=0;i<n_sub_list;i++) free(regions[i]);
         free(regions);
         return(siter);

      } else {  // no wildcard
         sequential_iter_t *siter = create_sequential_iter(t);
         add_to_sequential_iter ( siter, ti_querys_2d(t,reg) );
         return(siter);
      }
   }else {  // 1d query (let's regard it as 1d query on 1d index for now)
         sequential_iter_t *siter = create_sequential_iter(t);
         add_to_sequential_iter ( siter, ti_querys_2d(t,reg) );
         return(siter);
   }
}



ti_iter_t ti_queryi(pairix_t *t, int tid, int beg, int end)
{
        return ti_queryi_2d(t,tid,beg,end,-1,-1);
}

sequential_iter_t *ti_queryi_general(pairix_t *t, int tid, int beg, int end)
{
        return ti_queryi_2d_general(t,tid,beg,end,-1,-1);
}

ti_iter_t ti_queryi_2d(pairix_t *t, int tid, int beg, int end, int beg2, int end2)
{
	if (tid < 0) return ti_iter_first();
	if (ti_lazy_index_load(t) != 0) return 0;
	return ti_iter_query(t->idx, tid, beg, end, beg2, end2);
}

sequential_iter_t *ti_queryi_2d_general(pairix_t *t, int tid, int beg, int end, int beg2, int end2)
{
    sequential_iter_t *siter = create_sequential_iter(t);
    add_to_sequential_iter ( siter, ti_queryi_2d(t,tid,beg,end,beg2,end2) );
    return(siter);
}

ti_iter_t ti_querys(pairix_t *t, const char *reg)
{
        return ti_querys_2d(t,reg);
}

sequential_iter_t *ti_querys_general(pairix_t *t, const char *reg)
{
    sequential_iter_t *siter = create_sequential_iter(t);
    add_to_sequential_iter ( siter, ti_querys(t, reg) );
    return(siter);
}


ti_iter_t ti_querys_2d(pairix_t *t, const char *reg)
{
	int tid, beg, end, beg2, end2;
	if (reg == 0) return ti_iter_first();
	if (ti_lazy_index_load(t) != 0) return 0;
	if (ti_parse_region2d(t->idx, reg, &tid, &beg, &end, &beg2, &end2) < 0) return 0;
	return ti_iter_query(t->idx, tid, beg, end, beg2, end2);
}

sequential_iter_t *ti_querys_2d_multi(pairix_t *t, const char **regs, int nRegs)
{
    sequential_iter_t *siter = create_sequential_iter(t);
    int i;
    for(i=0;i<nRegs;i++){
      add_to_sequential_iter ( siter, ti_querys_2d(t, regs[i]) );
    }
    return(siter);
}

ti_iter_t ti_query(pairix_t *t, const char *name, int beg, int end)
{
	int tid;
	if (name == 0) return ti_iter_first();
	// then need to load the index
	if (ti_lazy_index_load(t) != 0) return 0;
	if ((tid = ti_get_tid(t->idx, name)) < 0) return 0;
	return ti_iter_query(t->idx, tid, beg, end, -1, -1);
}

sequential_iter_t *ti_query_general(pairix_t *t, const char *name, int beg, int end)
{
    sequential_iter_t *siter = create_sequential_iter(t);
    add_to_sequential_iter (siter, ti_query(t, name, beg, end));
    return(siter);
}

ti_iter_t ti_query_2d(pairix_t *t, const char *name, int beg, int end, const char *name2, int beg2, int end2)
{
	int tid;
        char namepair[1000], *str_ptr;
        char region_split_character = get_region_split_character(t);

        strcpy(namepair,name);
        str_ptr = namepair + strlen(namepair);
        *str_ptr = region_split_character;
        str_ptr++;
        strcpy(str_ptr,name2);

	if (name == 0) return ti_iter_first();
	// then need to load the index
	if (ti_lazy_index_load(t) != 0) return 0;
	if ((tid = ti_get_tid(t->idx, namepair)) < 0) return 0;
	return ti_iter_query(t->idx, tid, beg, end, beg2, end2);
}

sequential_iter_t *ti_query_2d_general(pairix_t *t, const char *name, int beg, int end, const char *name2, int beg2, int end2)
{
    sequential_iter_t *siter = create_sequential_iter(t);
    add_to_sequential_iter (siter, ti_query_2d(t, name, beg, end, name2, beg2, end2));
    return(siter);
}

int ti_querys_tid(pairix_t *t, const char *reg)
{
        return ti_querys_2d_tid(t,reg);
}

int ti_querys_2d_tid(pairix_t *t, const char *reg)
{
	int tid, beg, end, beg2, end2, parse_err;
	if (reg == 0) return -3;  // null region
	if (ti_lazy_index_load(t) != 0) return -3; // index not loaded
	parse_err = ti_parse_region2d(t->idx, reg, &tid, &beg, &end, &beg2, &end2);
        if(tid != -1 && tid != -3 && parse_err == -1) tid = -2;  // -2 is parsing error.
        return(tid);  // -1 means chromosome (pair) doesn't exist.
}

int ti_query_tid(pairix_t *t, const char *name, int beg, int end)
{
	if (name == 0) return -3 ;
	// then need to load the index
	if (ti_lazy_index_load(t) != 0) return -3;
	if (ti_get_tid(t->idx, name) <0) return -1;
        if (beg > end) return -2;
        else return 0;
}

int ti_query_2d_tid(pairix_t *t, const char *name, int beg, int end, const char *name2, int beg2, int end2)
{
        char namepair[1000], *str_ptr;
        char region_split_character = t->idx->conf.region_split_character;

        strcpy(namepair,name);
        str_ptr = namepair + strlen(namepair);
        *str_ptr = region_split_character;
        str_ptr++;
        strcpy(str_ptr,name2);

	if (name == 0) return (int)ti_iter_first();
	// then need to load the index
	if (ti_lazy_index_load(t) != 0) return 0;
        if(ti_get_tid(t->idx, namepair) <0) return -1;
        if (beg > end) return -2;
        if (beg2 > end2) return -2;
        else return 0;
}


const char *ti_read(pairix_t *t, ti_iter_t iter, int *len)
{
	return ti_iter_read(t->fp, iter, len, 0);
}


// create an empty sequential iterator
sequential_iter_t *create_sequential_iter(pairix_t *t)
{
  sequential_iter_t *siter = malloc(sizeof(sequential_iter_t));
  siter->t = t;
  siter->n = 0;
  siter->curr = 0;
  siter->iter = NULL;
  return(siter);
}

// destroy a sequential iterator
void destroy_sequential_iter(sequential_iter_t *siter)
{
  int i;
  if(siter){
    for(i=0;i<siter->n;i++) ti_iter_destroy(siter->iter[i]);
    free(siter->iter);
    free(siter);
  }
}

// add an iter element to a sequential iterator, the size is dynamically incremented.
void add_to_sequential_iter(sequential_iter_t *siter, ti_iter_t iter)
{
  siter->n++;
  siter->iter=realloc(siter->iter, siter->n * sizeof(ti_iter_t));
  siter->iter[siter->n-1] = iter;
}

// create an empty merged_iter_t struct with predefined length
merged_iter_t *create_merged_iter(int n)
{
   int i;
   merged_iter_t *miter = malloc(sizeof(merged_iter_t));
   if(miter){
     if( miter->iu = calloc(n,sizeof(iter_unit*))) {
       miter->n = n;
       miter->first=1;
       for(i=0;i<n;i++) miter->iu[i] = calloc(1,sizeof(iter_unit));
     } else { fprintf(stderr,"Cannot allocate memory for iter_unit array in merged_iter_t\n"); }
     return(miter);
   } else {
      fprintf(stderr,"Cannot allocate memory for merged_iter_t\n");
      return(NULL);
   }
}

void destroy_merged_iter(merged_iter_t *miter)
{
  int i;
  if(miter && miter->n>0 && miter->iu){
    for(i=0;i<miter->n;i++){
      ti_iter_destroy(miter->iu[i]->iter);
      if(miter->iu[i]->len) free(miter->iu[i]->len);
      if(miter->iu[i]) free(miter->iu[i]);
    }
    free(miter->iu);
    free(miter);
  }
}

// fill in an iter_unit struct given a pointer
void create_iter_unit(pairix_t *t, ti_iter_t iter, iter_unit *iu)
{
   if(iu){
     iu->t=t; iu->iter=iter; iu->len=malloc(sizeof(int)); iu->s=NULL;
   }
}


int compare_iter_unit (const void *a, const void *b)
{
  iter_unit *aa=*(iter_unit**)a;
  iter_unit *bb=*(iter_unit**)b;
  if (aa == NULL || aa->s == NULL) {
    if (bb == NULL) return(0);
    else if (bb->s == NULL) return(0);
    else return(1); // non-NULL iter should go before a NULL iter.
  } else if(bb == NULL || bb->s == NULL) {
    return(-1);
  } else {
    int res = aa->iter->intv.beg - bb->iter->intv.beg;  // sort first by beg
    if (res == 0 && aa->iter->intv.beg2 && bb->iter->intv.beg2) return ( aa->iter->intv.beg2 - bb->iter->intv.beg2 );  // sort second by beg2 (skip if beg2 doesn't exist - 1D case)
    else return (res);
  }
}


// read method for merged_iter
const char *merged_ti_read(merged_iter_t *miter, int *len)
{
    iter_unit *tmp_iu;
    int i,k;
    char *s;
    char seqonly=1; // seqonly=1 forces ti_iter_read to only check tid not positions. (faster, given the position comparison is done within this function)

    if(!miter) { fprintf(stderr,"Null merged_iter_t\n"); return(NULL); }
    if(miter->n<=0) { fprintf(stderr,"No iter_unit lement in merged_iter_t\n"); return(NULL); }
    iter_unit **miu = miter->iu;

    // initi al sorting of the iterators based on their first entry
    if (miter->first){
      for(i=0;i<miter->n;i++) {
        miu[i]->s = ti_iter_read(miu[i]->t->fp, miu[i]->iter, miu[i]->len, seqonly);
      } // get first entry for each iter
      qsort((void*)(miu), miter->n, sizeof(iter_unit*), compare_iter_unit);  // sort by the first entry. finished iters go to the end.
      miter->first=0;
    }
    else if(miu[0]->s==NULL) {
      miu[0]->s = ti_iter_read(miu[0]->t->fp, miu[0]->iter, miu[0]->len, seqonly); // get next entry for the flushed iter

      //qsort((void*)(miu), miter->n, sizeof(iter_unit*), compare_iter_unit); // sort again

      // put it at the right place in the sorted iu array
      k=0;
      while( k < miter->n-1 && compare_iter_unit((void*)(miu), (void*)(miu + k + 1 )) >0 ) k++;
      if (k>=1) {
        tmp_iu = miu[0];
        for(i=1;i<=k;i++) miu[i-1] = miu[i];
        miu[k]= tmp_iu;
      }
    }
    if(miu[0]->iter==NULL) return(NULL);

    // flush the lowest
    s=miu[0]->s;
    miu[0]->s=NULL;

    *len = *(miu[0]->len);

    return ( s );
}


const char *sequential_ti_read(sequential_iter_t *siter, int *len)
{
    if(!siter) { fprintf(stderr,"Null sequential_iter_t\n"); return(NULL); }
    if(siter->n<=0) { fprintf(stderr,"No iter_unit lement in merged_iter_t\n"); return(NULL); }

    char *s = ti_iter_read(siter->t->fp,siter->iter[siter->curr], len, 0);
    while(s==NULL && siter->curr < siter->n - 1) {
      siter->curr++;
      s = ti_iter_read(siter->t->fp,siter->iter[siter->curr], len, 0);
    }
    return s;
}




//compare two strings, but different from strcmp.
//for a pair of strings 'chr1|chr2' vs 'chr10|chr13', it compares chr1 vs chr10 first and then do chr2 vs chr13. This results in an ordering different from strcmp-based sort, because 'chr10' comes before 'chr1|' whereas 'chr1' comes before 'chr10'.
int strcmp2d(const void* a, const void* b)
{
    //char aa[strlen(*(char**)a)], bb[strlen(*(char**)b)];
    //strcpy(aa, *(char**)a);
    //strcpy(bb, *(char**)b);
    int res;
    char c,d;
    char *aa=*(char**)a;
    char *bb=*(char**)b;
    char *a2,*b2;
    char *a_split = strchr(aa, global_region_split_character);
    char *b_split = strchr(bb, global_region_split_character);
    if(a_split && b_split) {  // 2D name
      c = a_split[0]; d=b_split[0];
      a2 = a_split+1; b2 = b_split+1;
      a_split[0]=0;
      b_split[0]=0;
      if ((res = strcmp(aa, bb))==0) res =strcmp(a2, b2) ;
      a_split[0]=c; b_split[0]=d;
      return (res) ;
    } else if(!a_split && !b_split) {   // 1D name
       return strcmp(aa, bb);
    } else {  // weird mix?
       fprintf(stderr, "Warning: Mix of 1D and 2D indexed files? (%s vs %s)\n",aa,bb);
       return strcmp(aa, bb);
    }
}


// same as strcmp, argument types modified to be compatible with qsort
int strcmp1d(const void* a, const void* b)
{
   return( strcmp(*(char**)a, *(char**)b) );
}

pairix_t *load_from_file(char *fn)
{
      char *fnidx = calloc(strlen(fn) + 5, 1);
      strcat(strcpy(fnidx, fn), ".px2");
      pairix_t *tb = ti_open(fn, fnidx);
      free(fnidx);
      if(tb){
         tb->idx = ti_index_load(fn);
      } else return 0;
      return(tb);
}


// The returned array must be freed later.
char** get_unique_merged_seqname(pairix_t **tbs, int n, int *pn_uniq_seq)
{
     int i,j, k, len, n_seq_list=0;
     char **conc_seq_list=NULL;
     char **seqnames=NULL;
     if(n<=0) { fprintf(stderr, "Null pairix_t array\n"); return(0); }
     for(i=0;i<n;i++){
      if(tbs[i] && tbs[i]->idx){
            seqnames = (char**)ti_seqname(tbs[i]->idx,&len);
            if(seqnames){
              conc_seq_list = realloc(conc_seq_list, sizeof(char*)*(n_seq_list+len));
              for(k=0,j=n_seq_list;j<n_seq_list+len;j++,k++)
                conc_seq_list[j] = seqnames[k];
              n_seq_list += len;
              free(seqnames);
            }else { fprintf(stderr, "Cannot retrieve seqnames.\n"); return(0); }
      } else {
         for(j=0;j<i;j++) ti_close(tbs[j]);
         if(conc_seq_list) free(conc_seq_list);
         fprintf(stderr,"Not all files can be open.\n");
         return (0);
      }
    }

    if(conc_seq_list){
      // given an array, do sort|uniq, but doing it as if sorting by two chromosome columns (e.g. by chr1 first then chr2) rather than by a single merged chromosome pair string (e.g. 'chr1|chr2')
      qsort(conc_seq_list, n_seq_list, sizeof(char*), strcmp2d);  // This part does the sorting. see strcmp2d for more details.
      char **uniq_seq_list = uniq(conc_seq_list, n_seq_list, pn_uniq_seq);
      free(conc_seq_list);
      return ( uniq_seq_list );
    } else { fprintf(stderr,"Null concatenated seq list\n"); return(0); }
}


// given a chromosome for mate1 (seq1='chr1') return the array containing all seqpairs matching seq1 ('chr1', 'chr2', ... for seqpairs 'chr1|chr1', 'chr1|chr2', ... )
// the returned subarray contains copies of seq2 sequences (need to be freed later)
char **get_seq2_list_for_given_seq1(char *seq1, char **seqpair_list, int n_seqpair_list, int *pn_sub_list)
{
    int i,k;
    char *b_split, b;
    char **sublist;

    // first round, count the number
    k=0;
    for(i=0;i<n_seqpair_list;i++){
      b_split = strchr(seqpair_list[i], global_region_split_character);
      b = b_split[0];
      b_split[0] = 0;
      if ( strcmp(seqpair_list[i], seq1)==0 ) k++;
      b_split[0] = b;
    }
    *pn_sub_list = k;

    // second round, get the list of pointers
    sublist = malloc((*pn_sub_list)*sizeof(char*));
    k=0;
    for(i=0;i<n_seqpair_list;i++){
      b_split = strchr(seqpair_list[i], global_region_split_character);
      b = b_split[0];
      b_split[0] = 0;
      if ( strcmp(seqpair_list[i], seq1)==0 ) {
         sublist[k] = malloc((strlen(b_split+1)+1)*sizeof(char));
         strcpy(sublist[k], b_split+1);
         k++;
      }
      b_split[0] = b;
    }
    assert (k = *pn_sub_list);

    return(sublist);
}


// given a chromosome for mate2 (seq2='chr1') return the array containing all seqpairs matching seq2 ('chr1','chr2', ... for seqpairs 'chr1|chr1', 'chr2|chr1', ... )
// the returned subarray contains copies of seq1 sequences (need to be freed later)
char **get_seq1_list_for_given_seq2(char *seq2, char **seqpair_list, int n_seqpair_list, int *pn_sub_list)
{
    int i,k;
    char *b_split;
    char **sublist;

    // first round, count the number
    k=0;
    for(i=0;i<n_seqpair_list;i++){
      b_split = strchr(seqpair_list[i], global_region_split_character);
      if ( strcmp(b_split+1, seq2)==0 ) k++;
    }
    *pn_sub_list = k;

    // second round, get the list of pointers
    sublist = malloc((*pn_sub_list)*sizeof(char*));
    k=0;
    for(i=0;i<n_seqpair_list;i++){
      b_split = strchr(seqpair_list[i], global_region_split_character);
      if ( strcmp(b_split+1, seq2)==0 ) {
         *b_split=0;
         sublist[k] = malloc((strlen(seqpair_list[i])+1)*sizeof(char));
         strcpy(sublist[k], seqpair_list[i]);
         *b_split =  global_region_split_character;
         k++;
      }
    }
    assert (k = *pn_sub_list);

    return(sublist);
}

/* convert string 'region1|region2' to 'region2|region1' */
char *flip_region ( char* s, char region_split_character) {
    int l, i, l2, split_pos;
    l = strlen(s);
    char *s_flp = calloc(l+1, sizeof(char));
    for(i = 0; i != l; i++) if( s[i] == region_split_character) break;
    s[i]=0;
    split_pos = i;
    l2 = l-1-i;
    strcpy(s_flp, s + i + 1);
    s_flp[l2] = region_split_character;
    strcpy(s_flp + l2 + 1, s);
    s[i]= region_split_character;
    return(s_flp);
}

// given a chromosome for mate1 (seq1='chr1') return the array containing all seqpairs matching seq1 ('chr1|chr1', 'chr1|chr2', ... )
// the returned subarray contains pointers to the original seqpair_list elements.
char **get_sub_seq_list_for_given_seq1(char *seq1, char **seqpair_list, int n_seqpair_list, int *pn_sub_list)
{
    int i,k;
    char *b_split, b;
    char **sublist;

    // first round, count the number
    k=0;
    for(i=0;i<n_seqpair_list;i++){
      b_split = strchr(seqpair_list[i], global_region_split_character);
      b = b_split[0];
      b_split[0] = 0;
      if ( strcmp(seqpair_list[i], seq1)==0 ) k++;
      b_split[0] = b;
    }
    *pn_sub_list = k;

    // second round, get the list of pointers
    sublist = malloc((*pn_sub_list)*sizeof(char*));
    k=0;
    for(i=0;i<n_seqpair_list;i++){
      b_split = strchr(seqpair_list[i], global_region_split_character);
      b = b_split[0];
      b_split[0] = 0;
      if ( strcmp(seqpair_list[i], seq1)==0 ) { sublist[k] = seqpair_list[i]; k++; }
      b_split[0] = b;
    }
    assert (k = *pn_sub_list);

    return(sublist);
}


// given a chromosome for mate2 (seq2='chr1') return the array containing all seqpairs matching seq2 ('chr1|chr1', 'chr2|chr1', ... )
// the returned subarray contains pointers to the original seqpair_list elements.
char **get_sub_seq_list_for_given_seq2(char *seq2, char **seqpair_list, int n_seqpair_list, int *pn_sub_list)
{
    int i,k;
    char *b_split;
    char **sublist;

    // first round, count the number
    k=0;
    for(i=0;i<n_seqpair_list;i++){
      b_split = strchr(seqpair_list[i], global_region_split_character);
      if ( strcmp(b_split+1, seq2)==0 ) k++;
    }
    *pn_sub_list = k;

    // second round, get the list of pointers
    sublist = malloc((*pn_sub_list)*sizeof(char*));
    k=0;
    for(i=0;i<n_seqpair_list;i++){
      b_split = strchr(seqpair_list[i], global_region_split_character);
      if ( strcmp(b_split+1, seq2)==0 ) { sublist[k] = seqpair_list[i]; k++; }
    }
    assert (k = *pn_sub_list);

    return(sublist);
}


// given a list of seqpairs (either unique or non-unique), return an sorted array of unique seq1 list
// the returned array must be freed later. (first for the string elements and then itself)
char **get_seq1_list_from_seqpair_list(char** seqpair_list, int n_seqpair_list, int *pn_seq1)
{
    if(seqpair_list){
        char *b_split,b;
        char *seq1_list[n_seqpair_list];
        char **uniq_seq1_list=NULL;
        char *seqpair;
        int i;

        // extract seq1 from all seqpairs in the seqpair_list
        for(i=0;i<n_seqpair_list;i++){
          seqpair = seqpair_list[i];
          b_split = strchr(seqpair, global_region_split_character);
          b = b_split[0];
          b_split[0] = 0;
          seq1_list[i] = malloc((strlen(seqpair)+1)*sizeof(char));
          strcpy(seq1_list[i], seqpair);
          b_split[0] = b;
        }

        // uniquefy
        qsort(seq1_list, n_seqpair_list, sizeof(char*), strcmp1d); // sort first
        uniq_seq1_list = uniq(seq1_list, n_seqpair_list, pn_seq1); // uniq

        // free the non-unique seq1 list
        for(i=0;i<n_seqpair_list;i++) free(seq1_list[i]);

        // return the uniq list
        return(uniq_seq1_list);

    } else {
      fprintf(stderr, "Null seqpair list\n");
      return(0);
    }
}

char **uniq(char** seq_list, int n_seq_list, int *pn_uniq_seq)
{
    int k,i,prev_i;
    char **uniq_seq_list=NULL;

    // first round, just simulate uniquifying to get the length of the output array.
    k=0; prev_i=0;
    for(i=1;i<n_seq_list;i++){
      if ( strcmp(seq_list[i], seq_list[prev_i])==0) continue;   // This part does the uniquifying
      else {
        k++;
        prev_i=i;
      }
    }
    *pn_uniq_seq=k+1;
    fprintf(stderr,"(total %d unique seq names)\n",*pn_uniq_seq);

    // second round, allocate memory and actually create an array containing uniquified array
    if( (uniq_seq_list = malloc((*pn_uniq_seq)*sizeof(char*))) != NULL ) {
      k=0; prev_i=0;
      uniq_seq_list[0] = malloc((strlen(seq_list[0])+1)*sizeof(char));
      strcpy(uniq_seq_list[0],seq_list[0]);
      for(i=1;i<n_seq_list;i++){
        if ( strcmp(seq_list[i], seq_list[prev_i])==0) continue;  // uniquifying
        else {
          ++k;
          uniq_seq_list[k] = malloc((strlen(seq_list[i])+1)*sizeof(char));
          strcpy(uniq_seq_list[k],seq_list[i]);
          prev_i=i;
        }
      }
      assert(k+1==*pn_uniq_seq);
    } else { fprintf(stderr, "Cannot allocate memory for unique_seq_list\n"); return(0); }

    return(uniq_seq_list);
}


char get_region_split_character(pairix_t *t)
{
    return(t->idx->conf.region_split_character);
}


sequential_iter_t *querys_2D_wrapper(pairix_t *tb, const char *reg, int flip)
{
    int tid_test;
    sequential_iter_t *result;

    tid_test = ti_querys_tid(tb, reg);
    if (tid_test == -1) {
        char *reg2 = flip_region(reg, get_region_split_character(tb));
        int tid_test_rev = ti_querys_tid(tb, reg2);
        if (tid_test_rev != -1 && tid_test_rev != -2 && tid_test_rev != -3) {
            result = ti_querys_2d_general(tb, reg2);
            free(reg2);
            if (flip == 1){
                if (result == NULL) {
                   fprintf(stderr, "Cannot find matching chromosome pair. Check that chromosome naming conventions match between your query and input file.");
                   return(NULL);
                }else{
                    return(result);
                }
            }
            else{
                fprintf(stderr, "Cannot find matching chromosome pair. Check that chromosome naming conventions match between your query and input file. You may wish to also automatically test chromsomes in flipped order. To do this, include 1 as the last argument.");
                return(NULL);
            }
        }
        else free(reg2);
    }
    else if (tid_test == -2){
        fprintf(stderr, "The start coordinate must be less than the end coordinate.");
        return(NULL);
    }
    else if (tid_test == -3){
        fprintf(stderr, "The specific cause could not be found. Please adjust your arguments.");
        return(NULL);
    }

    result = ti_querys_2d_general(tb, reg);
    return(result);  // result may be null but that's okay
}

int get_nblocks(ti_index_t *idx, int tid, BGZF *fp)
{
    ti_iter_t iter = ti_iter_query(idx, tid, 0, 1<<MAX_CHR, 0, 1<<MAX_CHR);
    int64_t start_block_address = iter->off[0].u>>16;  // in bytes
    int64_t end_block_offset = iter->off[0].v;
    int nblocks=0;
    int64_t curr_off = start_block_address<<16;
    do {
      int block_length = bgzf_block_length(fp, curr_off);
      nblocks++;
      curr_off += block_length<<16;
    } while(curr_off <= end_block_offset);
    ti_iter_destroy(iter);

    return((int)nblocks);
}


// returns 1 if triangle
// returns 0 if not a triangle
// returns -1 if no chrom (pairs) is found in file
// returns -2 if the file is 1D-indexed (not applicable)
int check_triangle(ti_index_t *idx)
{
    if(ti_get_sc2(idx) == -1) return(-2);  // not a 2d file (not applicable)

    int len;
    const char **seqnames = ti_seqname(idx,&len);
    if(seqnames){
      int i;
      for(i=0;i<len;i++){
        const char *reg2 = flip_region(seqnames[i], ti_get_region_split_character(idx));
        if(strcmp(seqnames[i], reg2)!=0)
          if(ti_get_tid(idx, reg2)!=-1) { free(seqnames); free(reg2); return(0); }  // not a triangle
        free(reg2);
      }
      free(seqnames);
      return(1);  // is a triangle
    } else return(-1);  // no chromosome (pairs) found in file 

}


