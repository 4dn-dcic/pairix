#!/bin/bash

PATH=./bin:$PATH

if [ $VALGRIND_TEST_ON -eq 1 ]; then
  VALGRIND="valgrind --error-exitcode=42 --leak-check=full"
else
  VALGRIND=""
fi

## 2D
echo "test 1"
$VALGRIND pairix samples/merged_nodup.tab.chrblock_sorted.txt.gz '10:1-1000000|20' > log1
gunzip -c samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2=="10" && $3>=1 && $3<=1000000 && $6=="20"' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 1 failed"
  return 1;
fi

echo "test 1b"
$VALGRIND pairix -a samples/merged_nodup.tab.chrblock_sorted.txt.gz '10:1-1000000|20' > log1
gunzip -c samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2=="10" && $3>=1 && $3<=1000000 && $6=="20"' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 1b failed"
  return 1;
fi

echo "test 1c"
$VALGRIND pairix -a samples/merged_nodup.tab.chrblock_sorted.txt.gz '20|10:1-1000000' > log1
gunzip -c samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2=="10" && $3>=1 && $3<=1000000 && $6=="20"' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 1c failed"
  return 1;
fi

echo "test 1d"
$VALGRIND pairix samples/test_4dn.pairs.gz 'chr22:50000000-60000000' > log1
$VALGRIND pairix samples/test_4dn.pairs.gz 'chr22:50000000-60000000|chr22:50000000-60000000' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 1d failed"
  return 1;
fi

echo "test 1e"
$VALGRIND pairix samples/test_4dn.pairs.gz 'chrY:1-2000000' > log1
$VALGRIND pairix samples/test_4dn.pairs.gz 'chrY:1-2000000|chrY:1-2000000' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 1e failed"
  return 1;
fi

echo "test 2"
$VALGRIND pairix samples/merged_nodup.tab.chrblock_sorted.txt.gz '10:1-1000000|20:50000000-60000000' > log1
gunzip -c samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2=="10" && $3>=1 && $3<=1000000 && $6=="20" && $7>=50000000 && $7<=60000000' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 2 failed"
  return 1;
fi

echo "test 3"
$VALGRIND pairix samples/merged_nodup.tab.chrblock_sorted.txt.gz '1:1-10000000|20:50000000-60000000' '3:5000000-9000000|X:70000000-90000000' > log1
gunzip -c samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2=="1" && $3>=1 && $3<=10000000 && $6=="20" && $7>=50000000 && $7<=60000000' > log2
gunzip -c samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2=="3" && $3>=5000000 && $3<=9000000 && $6=="X" && $7>=70000000 && $7<=90000000' >> log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 3 failed"
  return 1;
fi

echo "test 4"
$VALGRIND pairix samples/merged_nodup.tab.chrblock_sorted.txt.gz '*|1:0-100000' > log1
gunzip -c samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$6=="1" && $7>=0 && $7<=100000' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 4 failed"
  return 1;
fi

echo "test 5"
$VALGRIND pairix samples/merged_nodup.tab.chrblock_sorted.txt.gz '1:0-100000|*' > log1
gunzip -c samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2=="1" && $3>=0 && $3<=100000' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 5 failed"
  return 1;
fi


## 1D
echo "test 6"
$VALGRIND pairix -s1 -b2 -e2 -f samples/SRR1171591.variants.snp.vqsr.p.vcf.gz
$VALGRIND pairix samples/SRR1171591.variants.snp.vqsr.p.vcf.gz chr10:1-4000000 > log1
gunzip -c samples/SRR1171591.variants.snp.vqsr.p.vcf.gz | awk '$1=="chr10" && $2>=1 && $2<=4000000' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 6 failed"
  return 1;
fi


## 2D, space-delimited
echo "test 7"
$VALGRIND pairix samples/merged_nodups.space.chrblock_sorted.subsample1.txt.gz '10:1-1000000|20' > log1
gunzip -c samples/merged_nodups.space.chrblock_sorted.subsample1.txt.gz | awk '$2=="10" && $3>=1 && $3<=1000000 && $6=="20"' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 7 failed"
  return 1;
fi


## preset for pairs.gz
echo "test 8"
$VALGRIND pairix samples/test_4dn.pairs.gz 'chr10|chr20' > log1
gunzip -c samples/test_4dn.pairs.gz | awk '$2=="chr10" && $4=="chr20"' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test 8 failed"
  return 1;
fi

## linecount
echo "test linecount"
$VALGRIND pairix -n samples/test_4dn.pairs.gz > log1
gunzip -c samples/test_4dn.pairs.gz |wc -l | sed "s/ //g" > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "linecount test failed"
  return 1;
fi


## bgzf block count (currently no auto test for the accuracy of the result)
echo "test bgzf block count"
$VALGRIND pairix -B samples/test_4dn.pairs.gz

## check triangle
echo "test check triangle"
$VALGRIND pairix -Y samples/4dn.bsorted.chr21_22_only.pairs.gz
$VALGRIND pairix -Y samples/4dn.bsorted.chr21_22_only.nontriangle.pairs.gz
res=$(pairix -Y samples/4dn.bsorted.chr21_22_only.nontriangle.pairs.gz)
if [ "$res" != "The file is not a triangle." ]; then
  echo "test check triangle failed"
  return 1;
fi

echo "test check triangle #2"
res=$(pairix -Y samples/4dn.bsorted.chr21_22_only.pairs.gz)
if [ "$res" != "The file is a triangle." ]; then
  echo "test check triangle #2 failed"
  return 1;
fi

# test large chromosome
echo "test large chr"
$VALGRIND pairix samples/mock.largechr.pairs.gz 'chr21:800000000-900000000|chr22' > log1
gunzip -c  samples/mock.largechr.pairs.gz | awk '$2=="chr21" && $3>800000000 && $3<900000000 && $4=="chr22"' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test large chromosome failed"
  return 1;
fi


# test large chromosome
echo "test large chr2"
$VALGRIND pairix samples/mock.largechr.pairs.gz 'chr22:800000000-997027270|chr22' > log1
gunzip -c samples/mock.largechr.pairs.gz | awk '$2=="chr22" && $3>=800000000 && $3<=997027270  && $4=="chr22"' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test large chromosome2 failed"
  return 1;
fi

# test large chromosome
echo "test large chr3"
$VALGRIND pairix samples/mock.largechr.pairs.gz 'chr22:1073741820-1073741824|chr22' > log1
gunzip -c samples/mock.largechr.pairs.gz | awk '$2=="chr22" && $3>=1073741820 && $3<=1073741824 && $4=="chr22"' > log2
if [ ! -z "$(diff log1 log2)" ]; then
  echo "test large chromosome3 failed"
  return 1;
fi


## pairs_merger
echo "test 13"
$VALGRIND pairs_merger samples/merged_nodups.space.chrblock_sorted.subsample2.txt.gz samples/merged_nodups.space.chrblock_sorted.subsample3.txt.gz | bgzip -c > out.gz
# compare with the approach of concatenating and resorting.
chmod +x test/inefficient-merger-for-testing
test/inefficient-merger-for-testing . out2 merged_nodups samples/merged_nodups.space.chrblock_sorted.subsample2.txt.gz samples/merged_nodups.space.chrblock_sorted.subsample3.txt.gz
gunzip -f out2.bsorted.pairs.gz 
gunzip -f out.gz
if [ ! -z "$(diff out out2.bsorted.pairs)" ]; then 
  echo "test 13 failed"
  return 1;
fi
rm -f out out2.bsorted.pairs out2.pairs out.gz.px2 out2.bsorted.pairs.gz.px2


## streamer_1d
echo "test 14"
$VALGRIND streamer_1d samples/merged_nodups.space.chrblock_sorted.subsample2.txt.gz | bgzip -c > out.1d.pairs.gz
gunzip -c samples/merged_nodups.space.chrblock_sorted.subsample2.txt.gz | sort -t' ' -k2,2 -k3,3g | bgzip -c > out2.1d.pairs.gz
gunzip -f out.1d.pairs.gz
gunzip -f out2.1d.pairs.gz
if [ ! -z "$(diff out.1d.pairs out2.1d.pairs)" ]; then
  echo "test 14 failed"
  return 1;
fi
rm -f out.1d.pairs out2.1d.pairs


## sort-pairs test
echo "test sort-pairs"
. util/sort-pairs.sh samples/tiny_unsorted.pairs.gz out
if [ ! -z "$(gunzip -c out.pairs.gz | diff - test/files/output_sort_tiny_unsorted.pairs)"]; then
  echo "test sort-pairs failed"
  return 1;
fi
rm -f out.pairs.gz*

