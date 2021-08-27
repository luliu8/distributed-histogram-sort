# README

Distributed Memory Histogram Sort

## What this is
Histogram sorting in a distributed memory environment using MPI parallelization.\
Benchmark performance on linux cluster using slurm. \
This is a project for the Parallel Programming course from UIUC

## Overall Algorithm 
1. balance data by redistributing it across processes without worrying about sorted-ness,
2. iteratively find a histogram with (roughly) equally counts in its bins by varying bin boundaries.
3. send data to the right ranks, so any data on rank i should be smaller than the smallest data item on rank i+1 (i.e. data should be globally sorted).
4. finally sort data in each rank locally. 

## implementation details and benchmark result
see: /writeup/proj_lul8.pdf