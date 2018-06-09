#!/bin/bash
./build.sh
echo "-------------------PREFIX-------------------"
./build/bin/main SHW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 100 1
echo
echo "-------------------GLOBAL-------------------"
./build/bin/main NW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 100 1
echo
echo "-------------------INFIX--------------------"
./build/bin/main HW ./test_data/test_pattern.fasta ./test_data/test_text.fasta 1  

dir=$(pwd)
mkdir -p results
cd /home/donik/Git/edlib/test_data
echo $dir

bash perf_tests2${1}.sh > "$dir/results/results${1}.txt" && perl $dir/resultsParser.pl $dir/results/results${1}.txt > "$dir/results/results${1}.csv"