#include <stdio.h>
#include <ctime>
#include <string>
#include <vector>
#include <unordered_map>
#include "FastaReader.hpp"
#include "Algorithms.hpp"

using namespace std;

int main (int argc, char** argv)
{
    vector<std::vector<char>> R;
    vector<std::vector<char>> B;

    readFastaSequences(argv[2], &R); readFastaSequences(argv[3], &B);
    vector<unsigned char> Rt(R.back().size()); vector<unsigned char> Bt(B.back().size());
    string alphabet = transformSequences(&(R.back())[0], (R.back()).size(), &(B.back())[0], (B.back()).size(), Rt, Bt);
    
    EqualityDefinition equality(alphabet);

    std::unordered_map<L, int, Hasher, EqualFn> D;
    int distance;
    bool global = string(argv[1]) == "NW";
    bool infix = string(argv[1]) == "HW";

    clock_t start = clock();

    for(int i=0; i<atoi(argv[4]); i++){
    if (global)
        distance = findAlgimentWithLowestKGLOBAL(Rt, Bt, equality);
    else if (infix)
        distance = findAlgimentWithLowestKINFIX(Rt, Bt, equality);
    else
        distance = findAlgimentWithLowestKPREFIX(Rt, Bt, equality);
    }
    printf("#0: %d\n", distance);

    clock_t finish = clock();
    double cpuTime = ((double)(finish-start))/CLOCKS_PER_SEC;
    printf("Cpu time of searching: %lf\n", cpuTime);

}