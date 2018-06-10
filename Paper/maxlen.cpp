// Example program
#include <iostream>
#include <string>

#include <vector>
#include <stdio.h>

void maxlength(const std::vector<unsigned int>& R, std::vector<unsigned int>& D, int m)
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

int main()
{
      int m = 12;
  std::vector<unsigned int> D(144);

  const unsigned int x[12] = {1,2,3,2,3,4,4,2,3,4,2,3};
    const std::vector<unsigned int> R({1,2,3,2,3,4,4,2,3,4,2,3});
    maxlength(R,D,m);
    for(int i=0;i<m;i++){
    for(int j=0;j<m;j++){
        if(j==m-1) printf("%d",D[i*m+j]);
        else printf("%d,",D[i*m+j]);
    }
    printf("\n");
}
}