#!/bin/bash

mkdir -p build
g++ -std=c++11 -o build/alg2 Algoritam2.cpp FastaReader.cpp EqualityDefinition.cpp 
g++ -std=c++11 -o build/alg3 Algoritam3.cpp Maxlength.cpp FastaReader.cpp EqualityDefinition.cpp 