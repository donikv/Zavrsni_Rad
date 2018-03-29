#include <vector>
#include <unordered_map>
#include <stdio.h>

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


bool operator==(const L& lhs ,const L& rhs)
{
        return lhs.d==rhs.d && lhs.e==rhs.e;
}

class Hasher
{
public:
  size_t operator() (L const& k) const
  {
    return (std::hash<std::string>()(std::to_string(k.d)+std::to_string(k.e)));
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

inline int max(int i1, int i2, int i3)
{
    //printf("%d,%d,%d\n", i1,i2,i3);
    if(i1>=i2 && i1>=i3) return i1;
    if(i2>=i1 && i2>=i3) return i2;
    return i3;
}

inline void printD(std::unordered_map<L, int, Hasher, EqualFn>& D)
{
    for(const auto l : D) printf("L(%d, %d) = %d\n", l.first.d, l.first.e, l.second);
}


bool algorithm(const std::vector<nucleotide>& R, const std::vector<nucleotide>& B, std::unordered_map<L, int, Hasher, EqualFn>& D, int m, int n, int bStart, int k)
{
    for (int d = -(k); d<=k; d++){
        D[L{d,abs(d)-2}] = -5;
        if(d<0) D[L{d, -d-1}] = -d-1;
        else D[L{d, d-1}] = -1;
    }

    unsigned int row = 0;
    for (int e = 0; e<=k; e++){
        for(int d = -e; d<=e; d++){
            printf("%d %d \n", d, e);
            if(d==-e){
                row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, -5);
            } else if (d==e) {
                row = max(D[L{d,e-1}]+1, -5, D[L{d-1,e-1}]);
            } else {
                row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, D[L{d-1,e-1}]);
            }
            //printf("while entrance: %d %d %d\n", row, d, row+d+bStart);
            
            while(R[row]==B[row+d+bStart] && row<m) row++;
            D[L{d,e}] = row;
            if(row == m){
                return true;
            }
        }
    }

    return false;
}

int main (void)
{
    std::vector<nucleotide> R = {A,G,C,G,C,T,T,G,C,T,G,C};
    std::vector<nucleotide> B = {A,G,T,C,G,C,C,G,C,T,G,C,T,G,C};

    int m = R.size();
    int n = B.size();
    printf("%d %d\n", m, n);

    std::unordered_map<L, int, Hasher, EqualFn> D;

    for(int i=0; i<n-m+3;i++){
        printf("i:%d %s\n", i, algorithm(R,B,D,m,n,i,3)==0 ? "NO" : "YES");
        printD(D);
        printf("\n");
    }
}
