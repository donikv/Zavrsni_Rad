#include <vector>
#include <stdio.h>

enum nucleotide : unsigned char {
    A = 0,
    C = 1,
    T = 2,
    G = 3
};

unsigned int min(unsigned int i1, unsigned int i2, unsigned int i3)
{
    if(i1<=i2 && i1<=i3) return i1;
    if(i2<=i1 && i2<=i3) return i2;
    return i3;
}

void algorithm(const std::vector<nucleotide>& R, const std::vector<nucleotide>& B, std::vector<unsigned int>& D, int m, int n)
{
    int row = n+1;
    D[0] = 0;
    for(int i=1;i<=n;i++) D[i] = i;
    for(int i=1;i<=m;i++) D[i*row] = i;

    for(int i=1;i<=m;i++){
        for(int j=1;j<=n;j++){
            D[i*row+j] = min(D[row*(i-1)+j]+1, D[row*i+j-1]+1, (R[i-1]==B[j-1]) ? D[row*(i-1)+j-1] : D[row*(i-1)+j-1]+1);
        }
    }
}

void printD(std::vector<unsigned int>& D, int m, int n)
{
    for(int i=0;i<=m;i++){
        for(int j=0;j<=n;j++){
            printf("%2d ", D[(n+1)*i+j]);
        }
        printf("\n");
    }
}

int main (void)
{
    std::vector<nucleotide> R = {A,G,C,G,C,T,T,G,C,T,G,C};
    std::vector<nucleotide> B = {A,G,T,C,G,C,C,G,C,T,G,C,T,G,C};

    int m = R.size();
    int n = B.size();
    printf("%d %d\n", m, n);

    std::vector<unsigned int> D((m+1)*(n+1));
    algorithm(R,B,D,m,n);
    printD(D, m, n);
}