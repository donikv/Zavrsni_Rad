#!/bin/bash

./build.sh || exit 1
echo "-------------------PREFIX-------------------"
./build/test SHW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1
echo
echo "-------------------GLOBAL-------------------"
./build/test HW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1
echo
echo "-------------------INFIX--------------------"
./build/test NW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1  