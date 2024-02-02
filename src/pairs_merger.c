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

int pairs_merger(char **fn, int n, BGZF *bzfp)  // pass bgfp if the result should be bgzipped. or pass NULL.
{
    pairix_t *tbs[n];
    int i,j;
    int reslen;
    int n_uniq_seq=0;
    char **uniq_seq_list=NULL;
    char *s=NULL;
    merged_iter_t *miter=NULL;
    ti_iter_t iter;

    // opening files and creating an array of pairix_t struct and prepare a concatenated seqname array
    fprintf(stderr,"Opening files...\n");
    for(i=0;i<n;i++)  {
       tbs[i] = load_from_file(fn[i]);
       if(i==0) global_region_split_character = get_region_split_character(tbs[i]);
       else if(global_region_split_character != get_region_split_character(tbs[i])){
           fprintf(stderr,"Merging is allowed only for files with the same region_split_character.\n"); return(1);
       }
    }

    // get a sorted unique seqname list
    fprintf(stderr,"creating a sorted unique seqname list...\n");
    uniq_seq_list = get_unique_merged_seqname(tbs, n, &n_uniq_seq);

    // loop over the seq_list (chrpair list) and merge
    if(uniq_seq_list){
      fprintf(stderr,"Merging...\n");
      for(i=0;i<n_uniq_seq;i++){
        miter = create_merged_iter(n);
        for(j=0;j<n;j++){
           iter = ti_querys_2d(tbs[j],uniq_seq_list[i]);
           create_iter_unit(tbs[j], iter, miter->iu[j]);
        }
        while ( ( s=merged_ti_read(miter,&reslen)) != NULL ) puts(s);
        destroy_merged_iter(miter); miter=NULL;
      }
      for(i=0;i<n;i++) ti_close(tbs[i]);
      for(i=0;i<n_uniq_seq;i++) free(uniq_seq_list[i]);
      free(uniq_seq_list);
      return(0);
    } else { fprintf(stderr,"Null unique seq list\n"); return(0); }
}

int main(int argc, char *argv[])
{
    int num_fn = argc -1;
    int i;
    //BGZF *bzfp;
    //int f_dst;

    if(argc==1){
       fprintf(stderr, "\n");
       fprintf(stderr, "Program: pairs_merger\n");
       fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
       //fprintf(stderr, "Usage:   pairs_merger <in1.pairs.gz> <in2.pairs.gz> <in3.pairs.gz> ... > out.pairs.gz\n\n");
       fprintf(stderr, "Usage:   pairs_merger <in1.pairs.gz> <in2.pairs.gz> <in3.pairs.gz> ... | bgzip -c  > out.pairs.gz\n\n");
       return(1);
    }

    char *fn_list[num_fn];
    for(i=0;i<num_fn;i++) {
       fn_list[i]=malloc(FILENAMEMAX*sizeof(char));
       strcpy(fn_list[i],argv[i+1]);
    }

    // write to stdout bgzip
    //f_dst = fileno(stdout);
    //bzfp = bgzf_dopen(f_dst, "w");

    // actually write merged pairs to bzfp stdout
    //int res = pairs_merger(fn_list, num_fn, bzfp);
    int res = pairs_merger(fn_list, num_fn, NULL);

    // close bgzf stream
    //if (bgzf_close(bzfp) < 0) fail(bzfp);
    //if (bgzf_close(bzfp) < 0){  fprintf(stderr,"Error: %d\n",bzfp->errcode); return(1); }


    for(i=0;i<num_fn;i++) {
       free(fn_list[i]);
    }
    return(res);
}
