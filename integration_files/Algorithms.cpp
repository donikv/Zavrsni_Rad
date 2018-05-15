#include "Algorithms.hpp"

int landauVishkinAlignAlgorithm(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, std::unordered_map<L, int, Hasher, EqualFn>& D, int m, int n, int bStart, int k, int nk, EqualityDefinition& equality, bool global, bool cigar, vector<char>& cigarVector)
{
    //bool standardCigar = false;
    unordered_map<L,vector<char>, Hasher, EqualFn> cigarDict;
    vector<char> cv;

    for (int d = -(k); d<=k; d++){
        if(d>=-nk && d<=nk && d!=0) continue;
        D[L{d,abs(d)-2}] = -5;
        if(d<0) D[L{d, -d-1}] = -d-1;
        else D[L{d, d-1}] = -1;

        if(cigar){
            cigarDict[L{d,abs(d)-1}] = cv;
            cigarDict[L{d,abs(d)-1}] = cv;
        }
    }

    unsigned int row = 0;
    int num = 0;
    for (int e = nk; e<=k; e++){
        for(int d = -e; d<=e; d++){
            if(d==-e){
                row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, -5, &num);
            } else if (d==e) {
                row = max(D[L{d,e-1}]+1, -5, D[L{d-1,e-1}], &num);
            } else {
                row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, D[L{d-1,e-1}], &num);
            }
            
            
            if(cigar){
                switch (num){
                    case 1:
                        cv = cigarDict[L{d,e-1}];
                        if(e!=0)
                            cv.push_back('I');
                        break;
                    case 2:
                        cv = cigarDict[L{d+1,e-1}];
                        cv.push_back('X');
                        break;
                    case 3:
                        cv = cigarDict[L{d-1,e-1}];
                        cv.push_back('D');
                        break;
                }
            }

            while(equality.areEqual(R[row], B[row+d+bStart]) && row<m) {
                if (cigar) cv.push_back('=');
                row++;
            }
            D[L{d,e}] = row;
            if (cigar) cigarDict[L{d,e}] = cv;

            if(row == m){
                if (cigar) cigarVector = cv;
                return global ? e+abs(n-m)+abs(d) : e;
            } else if (global && row+d == n) {
                return global ? e+abs(n-m)+abs(d) : e;
            }
        }
    }

    return -1;
}


int landauVishkinAlignPrefix(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar, string& cigarOutput)
{
    int m = R.size();
    int n = B.size();
    std::unordered_map<L, int, Hasher, EqualFn> D;
    int nk = 0;
    vector<char> cigarVector;
    int k = landauVishkinAlignAlgorithm(R,B,D,m,n,0,m, 0, equality, false, true, cigarVector);
    if (k!=-1) {
        if(cigar) {
            cigarOutput = standardCigarOutput(cigarVector);
        }  
        return k; 
    } 
    
    return -1;
}


string standardCigarOutput(vector<char>& cigarVector)
{
    ostringstream stream;
    stream << "POS: " << 0 << endl << "CIGAR: ";

    for(int i=0;i<cigarVector.size();i++) {
        char current = cigarVector[i];
        int j = 0;
        while(current==cigarVector[i] || (current=='=' && cigarVector[i]=='X') || (current=='X' && cigarVector[i]=='=')) {
            j++;
            i++;
        }
        i--;
        if(current=='='|| current=='X') {
            stream << j << "M";
        } else {
            stream << j << current;
        }
    }

    return stream.str();
}

string extendedCigarOutput(vector<char>& cigarVector)
{
    ostringstream stream;
    stream << "POS: " << 0 << endl << "CIGAR: ";

    for(int i=0;i<cigarVector.size();i++) {
        char current = cigarVector[i];
        int j = 0;
        while(current==cigarVector[i]) {
            j++;
            i++;
        }
        i--;
        stream << j << current;
    }

    return stream.str();
}