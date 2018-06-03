#include "Maxlength.hpp"

void maxlength(const std::vector<unsigned char>& R, std::vector<unsigned int>& D, int m, EqualityDefinition& equality)
{
    int row = m;

    for(int i=0;i<m;i++){
        for(int j=0;j<=i;j++){
            unsigned int f = 0;
            while(i+f<m && j+f<m && equality.areEqual(R[i+f], R[j+f])) f++;
            D[i*row+j] = f;
            D[row*j+i] = f;
        }
    }
}