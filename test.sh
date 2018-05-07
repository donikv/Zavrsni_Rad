#!/bin/bash
./build.sh
echo "-------------------PREFIX-------------------"
./build/bin/main SHW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1
echo
echo "-------------------GLOBAL-------------------"
./build/bin/main NW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1
echo
echo "-------------------INFIX--------------------"
./build/bin/main HW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1  