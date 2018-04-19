#!/bin/bash

./build.sh
echo "-------------------PREFIX-------------------"
./build/alg2 SHW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1
echo
echo "-------------------GLOBAL-------------------"
./build/alg2 HW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1
echo
echo "-------------------INFIX--------------------"
./build/alg3 ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1  