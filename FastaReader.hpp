#ifndef FASTA_READER_HPP
#define FASTA_READER_HPP

#include <stdio.h>
#include <vector>

using namespace std;

int readFastaSequences(const char* path, vector< vector<char> >* seqs);

#endif