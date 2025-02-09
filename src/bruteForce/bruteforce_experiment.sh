#!/bin/bash

# cd src

# CC bruteForce.cpp -o bf -fopenmp



export OMP_NUM_THREADS=1
./bf 2 8 >> ../../output/bruteforce/foutput.1.txt
./bf 2 9 >> ../../output/bruteforce/bfoutput.1.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.1.txt
./bf 2 11 >> ../../output/bruteforce/bfoutput.1.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.1.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.1.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.1.txt

export OMP_NUM_THREADS=2
./bf 2 8 >> ../../output/bruteforce/bfoutput.2.txt
./bf 2 9 >> ../../output/bruteforce/bfoutput.2.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.2.txt
./bf 2 11 >> ../../output/bruteforce/bfoutput.2.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.2.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.2.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.2.txt

export OMP_NUM_THREADS=4
./bf 2 8 >> ../../output/bruteforce/bfoutput.4.txt
./bf 2 9 >> ../../output/bruteforce/bfoutput.4.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.4.txt
./bf 2 11 >> ../../output/bruteforce/bfoutput.4.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.4.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.4.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.4.txt

export OMP_NUM_THREADS=8
./bf 2 8 >> ../../output/bruteforce/bfoutput.8.txt
./bf 2 9 >> ../../output/bruteforce/bfoutput.8.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.8.txt
./bf 2 11 >> ../../output/bruteforce/bfoutput.8.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.8.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.8.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.8.txt

export OMP_NUM_THREADS=16
./bf 2 8 >> ../../output/bruteforce/bfoutput.16.txt
./bf 2 9 >> ../../output/bruteforce/bfoutput.16.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.16.txt
./bf 2 11 >> ../../output/bruteforce/bfoutput.16.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.16.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.16.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.16.txt

export OMP_NUM_THREADS=32
./bf 2 8 >> .../../output/bruteforce/bfoutput.32.txt
./bf 2 9 >> ../../output/bruteforce/bfoutput.32.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.32.txt
./bf 2 11 >> .../../output/bruteforce/bfoutput.32.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.32.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.32.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.32.txt

export OMP_NUM_THREADS=64
./bf 2 8 >>../../output/bruteforce/bfoutput.64.txt
./bf 2 9 >> ../../output/bruteforce/bfoutput.64.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.64.txt
./bf 2 11 >> ../../output/bruteforce/bfoutput.64.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.64.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.64.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.64.txt


export OMP_NUM_THREADS=128
./bf 2 8 >> ../../output/bruteforce/bfoutput.128.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.128.txt
./bf 2 11 >> ../../output/bruteforce/bfoutput.128.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.128.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.128.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.128.txt



export OMP_NUM_THREADS=256
./bf 2 8 >> ../../output/bruteforce/bfoutput.256.txt
./bf 2 9 >> ../../output/bruteforce/bfoutput.256.txt
./bf 2 10 >> ../../output/bruteforce/bfoutput.256.txt
./bf 2 11 >> ../../output/bruteforce/bfoutput.256.txt
./bf 2 12 >> ../../output/bruteforce/bfoutput.256.txt
./bf 2 13 >> ../../output/bruteforce/bfoutput.256.txt
./bf 2 14 >> ../../output/bruteforce/bfoutput.256.txt









