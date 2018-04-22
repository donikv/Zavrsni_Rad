#!/bin/bash

mkdir -p build
g++ -std=c++11 -o build/test main.cpp Algorithms.cpp Maxlength.cpp FastaReader.cpp EqualityDefinition.cpp 