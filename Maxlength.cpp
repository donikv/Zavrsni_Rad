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

void printD(std::vector<unsigned int>& D, int m)
{
    for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
            printf("%2d ", D[(m)*i+j]);
        }
        printf("\n");
    }
}

int main (void)
{
    std::vector<char> R = {A,G,C,G,C,T,T,G,C,T,G,C};

    int m = R.size();
    printf("%d\n", m);

    std::vector<unsigned int> D((m)*(m));
    maxlength(R,D,m);
    printD(D, m);
}