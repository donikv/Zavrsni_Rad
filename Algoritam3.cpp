#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <string>
#include "FastaReader.cpp"

enum nucleotide : unsigned char {
    A = 0,
    C = 1,
    T = 2,
    G = 3
};

struct L {
    int d;
    int e;
};

struct Triple {
    unsigned int p;
    unsigned int c;
    unsigned int f;
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


int algorithm(const std::vector<char>& R, const std::vector<char>& B, const std::vector<unsigned int>& MAXLENGTH,int m, int n, int kmax, std::unordered_map<int, bool>& ret)
{
    for(int k = 0; k<kmax; k++)
    {
    int j = 0;
    std::vector<Triple> Sij;

    for(int i=0; i<n-m+k;i++){
        std::unordered_map<L, std::vector<Triple>, Hasher, EqualFn> lSeqMap;
        std::unordered_map<L, int, Hasher, EqualFn> D;

        for (int d = -(k); d<=k; d++){
            D[L{d,abs(d)-2}] = -5;
            if(d<0) D[L{d, -d-1}] = -d-1;
            else D[L{d, d-1}] = -1;
            lSeqMap[L{d,abs(d)-2}] = std::vector<Triple>();
            lSeqMap[L{d,abs(d)-1}] = std::vector<Triple>();
        }

        unsigned int row = 0;
        unsigned int l1;
        int num = 0;
        int e = 0;
        int d = -e;

        for (; e<=k; e++){
            d=-e;
            for(; d<=e; d++){
                if(d==-e){
                    row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, -5, &num);
                } else if (d==e) {
                    row = max(D[L{d,e-1}]+1, -5, D[L{d-1,e-1}], &num);
                } else {
                    row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, D[L{d-1,e-1}], &num);
                }

                l1 = row;
                while(row+d+i<=j){
                    unsigned int c=0;
                    unsigned int f=0;
                    for(const auto& t: Sij){
                        if(t.p+t.f>row+d+i && t.p==row+d+i){
                            f=t.f;
                            c=t.c;
                            break;
                        }
                    }
                    //printf("After for: c,f,row,maxlen %d,%d,%d,%d\n", c, f, row, MAXLENGTH[c*m+row]);

                    if(f>=1) {
                        if(f!=MAXLENGTH[c*m+row]) {
                            row += std::min(f,MAXLENGTH[c*m+row]);
                            goto inst5;
                        } else {
                            row += f;
                        }
                    } else {
                        if(R[row]!=B[row+d+i]){
                            goto inst5;
                        } else {
                            row++;
                        }
                    }
                }

                while(R[row]==B[row+d+i] && row<m) row++;

            inst5:
                D[L{d,e}] = row;

                L pickedL; //get L that was picked as the L that gives the maximum row
                if(num==1) { pickedL.d = d; pickedL.e = e-1; }
                else if(num == 2){ pickedL.d = d+1; pickedL.e = e-1; }
                else { pickedL.d = d-1; pickedL.e = e-1; }
                //printf("inst5 L(%d %d) = %d, Pl(%d %d), l1 = %d ",d, e, row, pickedL.d, pickedL.e, l1);

                std::vector<Triple> sequenceForCurrentL = lSeqMap[pickedL];
                //for(const auto& seq : sequenceForCurrentL) printf("(%d %d %d)", seq.p, seq.c, seq.f);
                //printf(" ");
                if((pickedL.d == d-1 || pickedL.d == d) && l1+d>0) sequenceForCurrentL.push_back(Triple{i+l1+d-1,0,0});
                if(row>l1) sequenceForCurrentL.push_back(Triple{i+l1+d, l1, row-l1});
                lSeqMap[L{d,e}] = sequenceForCurrentL;
                //for(const auto& seq : sequenceForCurrentL) printf("(%d %d %d)", seq.p, seq.c, seq.f);
                //printf("\n");

                if(row == m){
                    //ret[i] = true;
                    goto inst7;
                }
            }
        }

    inst7:
    //printD(D);
    //printf("\n");
    //printf("inst7 %d %d\n", row+d+i, j);
        //printf("2. p,c,f %d,%d,%d\n", i+l1+d-1,0,0);
        //printf("3. p,c,f %d,%d,%d\n", i+l1+d, l1, row-l1);
        if(i+row+d<=j) continue;
        j = i+row+d;

        for(int l = -k; l<=k; l++){
            std::vector<Triple> sequenceForCurrentL = lSeqMap[L{l,k}];
            if(sequenceForCurrentL.size()<=0) continue;

            if((sequenceForCurrentL.back().p+sequenceForCurrentL.back().f)>=j) {
                // j = sequenceForCurrentL.back().p+sequenceForCurrentL.back().f;
                Sij = sequenceForCurrentL;
                break;
            }
        }

        for(const auto& seq : Sij) 
        {
            printf("(%d %d %d)", seq.p, seq.c, seq.f);
        }
        printf("\n");
        int numInsertions = 0, numDeletions = 0, numChange = 0;
        for(int seqNum = 0; seqNum<Sij.size(); seqNum++) 
        {
            Triple seq = Sij[seqNum];
            if(!seqNum) { printf("POS: %d\nCIGAR: %dM", seq.p, seq.f); }
            else if(seq.f==0 && seq.c==0)
            {
                numInsertions++;
            }
            else
            {
                Triple previous = Sij[seqNum-1-numInsertions];
                //if (numInsertions!=0 && previous.c+previous.f+numInsertions != seq.c) printf("%dI", numInsertions);
                if (numInsertions!=0) printf("%dI", numInsertions);
                
                if(previous.c+previous.f == seq.c) printf("%dM", seq.f);
                else if(previous.c+previous.f+numInsertions == seq.c && previous.p+previous.f+numInsertions == seq.p) printf("%dM", seq.c-previous.c-previous.f);
                else printf("%dD", seq.c-previous.c-previous.f);
                numInsertions = 0;
            }
            //printf("(%d %d %d)", seq.p, seq.c, seq.f);
            
        }
        printf("\n");
        if(row == m) { printf("\n"); return k; } 
    }
    //for(const auto& i : ret) printf("%d %d\n", i.first,i.second);
    }
    return -1;
}

void maxlength(const std::vector<char>& R, std::vector<unsigned int>& D, int m)
{
    int row = m;

    for(int i=0;i<m;i++){
        for(int j=0;j<=i;j++){
            unsigned int f = 0;
            while(i+f<m && j+f<m && R[i+f]==R[j+f]) f++;
            D[i*row+j] = f;
            D[row*j+i] = f;
        }
    }
}

int main (int argc, char** argv)
{
    vector<std::vector<char>> R;
    vector<std::vector<char>> B;

    readFastaSequences(argv[1], &R); readFastaSequences(argv[2], &B); 

    int m = R.back().size();
    int n = B.back().size();
    //printf("%d %d\n", m, n);

    std::vector<unsigned int> MAXLENGTH((m)*(m));
    maxlength(R.back(),MAXLENGTH,m);

    std::unordered_map<int, bool> ret;

    clock_t start = clock();
    for(int i=1; i<atoi(argv[3]); i++) algorithm(R.back(),B.back(),MAXLENGTH,m,n,m,ret);
    printf("#0: %d\n", algorithm(R.back(),B.back(),MAXLENGTH,m,n,m,ret));
    //for(const auto& i : ret) printf("%d %s\n", i.first,i.second ? "YES" : "NO");
    //printf("\n");

    clock_t finish = clock();
    double cpuTime = ((double)(finish-start))/CLOCKS_PER_SEC;
    printf("Cpu time of searching: %lf\n", cpuTime);
}