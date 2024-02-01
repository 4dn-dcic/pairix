#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#include "bgzf.h"
#include "pairix.h"
#include "knetfile.h"

#define FILENAMEMAX 2000


// Uc->Up converter - convert a single 2D-sorted file into a 1D-sorted stream.
int stream_1d(char *fn)
{
    pairix_t *tb, **tbs_copies=NULL;
    int n_chrpairs, n_chr1, n_chr1pairs;
    char **chrpair_list, **chr1_list, **chr1pair_list;
    int i,j;
    char *s=NULL;
    merged_iter_t *miter=NULL;
    ti_iter_t iter;
    int reslen;

    tb = load_from_file(fn);
    if(tb==NULL) { fprintf(stderr,"file load failed\n"); return(1); }

    global_region_split_character = get_region_split_character(tb);

    chrpair_list = ti_seqname(tb->idx, &n_chrpairs);
    if(chrpair_list==NULL) { fprintf(stderr, "Cannot retrieve key list\n"); return(1); }
    chr1_list = get_seq1_list_from_seqpair_list(chrpair_list, n_chrpairs, &n_chr1);  // 'chr1','chr2',...
    if(chr1_list==NULL) { fprintf(stderr, "Cannot retrieve list of first chromosomes\n"); return(1); }

    for(i=0;i<n_chr1;i++){
       chr1pair_list = get_sub_seq_list_for_given_seq1(chr1_list[i], chrpair_list, n_chrpairs, &n_chr1pairs); // 'chr2|chr2', 'chr2|chr3' ... given chr2, this one is not necessarily a sorted list but it doesn't matter.
       miter = create_merged_iter(n_chr1pairs);
       tbs_copies= malloc(n_chr1pairs*sizeof(pairix_t*));
       for(j=0;j<n_chr1pairs;j++){
           tbs_copies[j] = load_from_file(fn);
           iter = ti_querys_2d(tbs_copies[j],chr1pair_list[j]);
           create_iter_unit(tbs_copies[j], iter, miter->iu[j]);
       }
       while ( (s=merged_ti_read(miter,&reslen)) != NULL ) puts(s);
       destroy_merged_iter(miter); miter=NULL;
       for(j=0;j<n_chr1pairs;j++) ti_close(tbs_copies[j]);
       free(tbs_copies); tbs_copies=NULL;
       free(chr1pair_list);
    }

    ti_close(tb);
    for(i=0;i<n_chr1;i++) free(chr1_list[i]);
    free(chr1_list);
    free(chrpair_list);

    return (0);
}

int main(int argc, char *argv[])
{
    //BGZF *bzfp;
    //int f_dst;

    if(argc==1){
       fprintf(stderr, "\n");
       fprintf(stderr, "Program: streamer_1d\n");
       fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
       fprintf(stderr, "Resorter (convert a file sorted by chr1-chr2-pos1-pos2 to a stream sorted by chr1-pos1)\n\n");
       fprintf(stderr, "Usage:   streamer_1d in.2d.pairs.gz > out.1d.pairs\n");
       fprintf(stderr, "Usage:   streamer_1d in.2d.pairs.gz | bgzip -c  > out.1d.pairs.gz\n\n");
       return(1);
    }

    char *fn;
    fn=malloc(FILENAMEMAX*sizeof(char));
    strcpy(fn,argv[1]);

    // write to stdout bgzip (This is slower, so don't use it)
    //f_dst = fileno(stdout);
    //bzfp = bgzf_dopen(f_dst, "w");

    // actually write merged pairs to bzfp stdout
    //int res = stream_1d(fn, bzfp);
    int res = stream_1d(fn);

    // close bgzf stream
    //if (bgzf_close(bzfp) < 0) fail(bzfp);
    //if (bgzf_close(bzfp) < 0){  fprintf(stderr,"Error: %d\n",bzfp->errcode); return(1); }

    free(fn);
    return(res);
}

