#include <vector>
#include <unordered_map>
#include <stdio.h>
#include "FastaReader.hpp"
#include "EqualityDefinition.hpp"


# define my_sizeof(type) ((char *)(&type+1)-(char*)(&type))

using namespace std;

enum nucleotide : unsigned char {
    A = 'A',
    C = 'C',
    T ='T',
    G ='G'
};

struct L {
    int d;
    int e;
};

struct Cigar {
    char type;
    unsigned int num;
};

bool operator==(const L& lhs ,const L& rhs)
{
        return lhs.d==rhs.d && lhs.e==rhs.e;
}

class Hasher
{
public:
  size_t operator() (L const& k) const
  {
    //return (std::hash<std::string>()(std::to_string(k.d)+std::to_string(k.e)));
    return k.d*100+k.e;
  }
};
class EqualFn
{
public:
  bool operator() (L const& t1, L const& t2) const
  {
    return t1==t2;
  }
};

inline int max(int i1, int i2, int i3, int* num)
{
    //printf("%d,%d,%d\n", i1,i2,i3);
    if(i1>=i2 && i1>=i3) {*num = 1; return i1;}
    if(i2>=i1 && i2>=i3) {*num = 2; return i2;}
    *num = 3;
    return i3;
}

inline void printD(std::unordered_map<L, int, Hasher, EqualFn>& D)
{
    for(const auto l : D) printf("L(%d, %d) = %d\n", l.first.d, l.first.e, l.second);
}


int algorithm(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, std::unordered_map<L, int, Hasher, EqualFn>& D, int m, int n, int bStart, int k, int nk, EqualityDefinition& equality)
{
    //bool standardCigar = false;  

    for (int d = -(k); d<=k; d++){
        if(d>=-nk && d<=nk) continue;
        D[L{d,abs(d)-2}] = -5;
        if(d<0) D[L{d, -d-1}] = -d-1;
        else D[L{d, d-1}] = -1;
    }

    unsigned int row = 0;
    int num = 0;
    for (int e = nk; e<=k; e++){
        //std::vector<Cigar> cigarOutput;
        for(int d = -e; d<=e; d++){
            //printf("%d %d \n", d, e);
            if(d==-e){
                row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, -5, &num);
            } else if (d==e) {
                row = max(D[L{d,e-1}]+1, -5, D[L{d-1,e-1}], &num);
            } else {
                row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, D[L{d-1,e-1}], &num);
            }
            //printf("while entrance: %d %d %d\n", row, d, row+d+bStart);
            // if (num==1) {
            //     if(standardCigar){
            //         if(cigarOutput.size() && cigarOutput.back().type == 'M') cigarOutput.back().num++;
            //         else cigarOutput.push_back(Cigar{'M', 1});
            //     } else {
            //         if(cigarOutput.size() && cigarOutput.back().type == 'X') cigarOutput.back().num++;
            //         else cigarOutput.push_back(Cigar{'X', 1});
            //     }
            // } else if (num==2) {
            //     if(cigarOutput.size() && cigarOutput.back().type == 'D') cigarOutput.back().num++;
            //     else cigarOutput.push_back(Cigar{'D', 1});
            // } else {
            //     if(cigarOutput.size() && cigarOutput.back().type == 'I') cigarOutput.back().num++;
            //     else cigarOutput.push_back(Cigar{'I', 1});
            // }
            // unsigned int prevRow = row;

            while(equality.areEqual(R[row], B[row+d+bStart]) && row<m) row++;
            D[L{d,e}] = row;

            // if(row-prevRow!=0) { 
            //     if(cigarOutput.size() && cigarOutput.back().type == 'M') cigarOutput.back().num+=(row-prevRow);
            //     else cigarOutput.push_back(Cigar{'M', row-prevRow});
            // } 

            if(row == m){
                // printf("POS: 0\nCIGAR: ");
                // for (auto& c : cigarOutput) printf("%d%c", c.num, c.type);
                printf("\n");
                return e;
            }
        }
    }

    return -1;
}

//prefix
int findAlgimentWithLowestKPREFIX(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality)
{
    int m = R.size();
    int n = B.size();
    std::unordered_map<L, int, Hasher, EqualFn> D;
    int nk = 0;
    for(int i=4;i>0;i/=2)
    {  
        if(nk != 0) nk = m/nk;
        int k = algorithm(R,B,D,m,n,0,m/i, nk, equality);
        if (k!=-1) return k;
        nk = i;
    }
    int k = algorithm(R,B,D,m,n,0,m/4, 0, equality);
    if (k!=-1) return k;
    return algorithm(R,B,D,m,n,0,m, m/2, equality);
}

int findAlgimentWithLowestKGLOBAL(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality)
{
    int m = R.size();
    int n = B.size();

    std::unordered_map<L, int, Hasher, EqualFn> D;

    return (findAlgimentWithLowestKPREFIX(R, B, equality) + n - m);
}

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
    bool global = string(argv[1]) == "HW";

    clock_t start = clock();

    for(int i=0; i<atoi(argv[4]); i++)
    if (global)
        distance = findAlgimentWithLowestKGLOBAL(Rt, Bt, equality);
    else
        distance = findAlgimentWithLowestKPREFIX(Rt, Bt, equality);
    printf("#0: %d\n", distance);

    clock_t finish = clock();
    double cpuTime = ((double)(finish-start))/CLOCKS_PER_SEC;
    printf("Cpu time of searching: %lf\n", cpuTime);

}
