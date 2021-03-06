#include "Algorithms.hpp"

int algorithm2(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, std::unordered_map<L, int, Hasher, EqualFn>& D, int m, int n, int bStart, int k, int nk, EqualityDefinition& equality, bool global, bool cigar, vector<char>& cigarVector)
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
            cigarDict[L{d,abs(d)-2}] = cv;
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

int algorithm2Global(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, std::unordered_map<L, int, Hasher, EqualFn>& D, int m, int n, int bStart, int k, int nk, EqualityDefinition& equality, bool cigar, vector<char>& cigarVector)
{
    //bool standardCigar = false;
    unordered_map<L,vector<char>, Hasher, EqualFn> cigarDict;
    vector<char> cv;
    for(int i = 0; i<bStart;i++) cv.push_back('I');

    for (int d = -(k); d<=k; d++){
        if(d>=-nk && d<=nk && d!=0) continue;
        D[L{d,abs(d)-2}] = -5;
        if(d<0) D[L{d, -d-1}] = -d-1;
        else D[L{d, d-1}] = -1;

        if(cigar){
            cigarDict[L{d,abs(d)-2}] = cv;
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

            while(equality.areEqual(R[row], B[row+d+bStart]) && row<m && row+d+bStart < n) {
                if (cigar) cv.push_back('=');
                row++;
            }
            D[L{d,e}] = row;
            if (cigar) cigarDict[L{d,e}] = cv;

            if(row == m){
                if (cigar) cigarVector = cv;
                return e+(n)-m-d + bStart;
            } else if (row+d+bStart == n) {
                return e-(n)+m+d + bStart;
            }
        }
    }

    return -1;
}

int algorithm3(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, const std::vector<unsigned int>& MAXLENGTH,int m, int n, int k, EqualityDefinition& equality, bool prefix, bool cigar, vector<char>& cigarVector)
{
    int j = 0;
    int e = 0;
    std::vector<Triple> Sij;
    unordered_map<L,vector<char>, Hasher, EqualFn> cigarDict;

    for(int i=0; i<(prefix ? 1 : n-m+k);i++){
        std::unordered_map<L, std::vector<Triple>, Hasher, EqualFn> lSeqMap;
        std::unordered_map<L, int, Hasher, EqualFn> D;
        vector<char> cv;

        for (int d = -(k); d<=k; d++){
            D[L{d,abs(d)-2}] = -5;
            if(d<0) D[L{d, -d-1}] = -d-1;
            else D[L{d, d-1}] = -1;
            lSeqMap[L{d,abs(d)-2}] = std::vector<Triple>();
            lSeqMap[L{d,abs(d)-1}] = std::vector<Triple>();

            if(cigar){
                cigarDict[L{d,abs(d)-2}] = cv;
                cigarDict[L{d,abs(d)-1}] = cv;
            }
        }

        unsigned int row = 0;
        unsigned int l1;
        int num = 0;
        e = 0;
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

                    if(f>=1) {
                        if(f != MAXLENGTH[c*m+row]) {
                            row += std::min(f,MAXLENGTH[c*m+row]);
                            if (cigar) {
                                for(int a=0;a<std::min(f,MAXLENGTH[c*m+row]);a++) cv.push_back('=');
                            }
                            goto inst5;
                        } else {
                            if (cigar) {
                                for(int a=0;a<f;a++) cv.push_back('=');
                            }
                            row += f;
                        }
                    } else {
                        if(!(equality.areEqual(R[row], B[row+d+i]))){
                            goto inst5;
                        } else {
                            if (cigar) cv.push_back('=');
                            row++;
                        }
                    }
                }

                while(equality.areEqual(R[row], B[row+d+i]) && row<m) { 
                    row++; 
                    if (cigar) cv.push_back('='); 
                } 

            inst5:
                D[L{d,e}] = row;
                if (cigar) cigarDict[L{d,e}] = cv;

                L pickedL; //get L that was picked as the L that gives the maximum row
                if(num==1) { pickedL.d = d; pickedL.e = e-1; }
                else if(num == 2){ pickedL.d = d+1; pickedL.e = e-1; }
                else { pickedL.d = d-1; pickedL.e = e-1; }

                std::vector<Triple> sequenceForCurrentL = lSeqMap[pickedL];

                if((pickedL.d == d-1 || pickedL.d == d) && l1+d>0) sequenceForCurrentL.push_back(Triple{i+l1+d-1,0,0});
                if(row>l1) sequenceForCurrentL.push_back(Triple{i+l1+d, l1, row-l1});
                lSeqMap[L{d,e}] = sequenceForCurrentL;

                if(row == m){
                    //ret[i] = true;
                    goto inst7;
                }
            }
        }

    inst7:
        if(i+row+d<=j) continue;
        j = i+row+d;

        L current_L;

        for(int l = -k; l<=k; l++){//TODO SMISLITI BOLJI NAČIN KOJI NIJE O(k)
            current_L = L{l,k};
            std::vector<Triple> sequenceForCurrentL = lSeqMap[current_L];
            if(sequenceForCurrentL.size()<=0) continue;

            if((sequenceForCurrentL.back().p+sequenceForCurrentL.back().f)>=j) {
                Sij = sequenceForCurrentL;
                break;
            }
        }
        if(row == m) { 
            if (cigar) cigarVector = cigarDict[current_L]; 
            return e; 
        }
    }
    return -1;
}

//prefix
int findAlgimentWithLowestKPREFIX(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar, string& cigarOutput)
{
    int m = R.size();
    int n = B.size();
    std::unordered_map<L, int, Hasher, EqualFn> D;
    int nk = 0;
    vector<char> cigarVector;

    for(int i=4;i>0;i/=2)
    {  
        if(nk != 0) nk = m/nk;
        int k = algorithm2(R,B,D,m,n,0,m, 0, equality, false, true, cigarVector);
        if (k!=-1) {
            if(cigar) {
                cigarOutput = standardCigarOutput(cigarVector);
            }  
            return k; 
        } 
        nk = i;
    }
    
    return -1;
}

int findAlgimentWithLowestKGLOBAL(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar, string& cigarOutput)
{
    int m = R.size();
    int n = B.size();
    std::unordered_map<L, int, Hasher, EqualFn> D;
    int nk = 0;
    vector<char> cigarVector;
    int mink = -1;

    for(int i=0;i<n;i++)
    {  
        if(i>mink && mink!=-1) break;
        int k = algorithm2Global(R,B,D,m,n,i,m, nk, equality, true, cigarVector);
        // if(cigar) {
        //     cigarOutput = standardCigarOutput(cigarVector);
        // }
        // return k;
        if (k!=-1) { 
            if(k<mink || mink == -1) {
                mink = k;
                if(cigar) {
                    cigarOutput = standardCigarOutput(cigarVector);
                } 
            }
        } 
    }
    return mink;
}

int findAlgimentWithLowestKINFIX(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar, string& cigarOutput)
{
    int m = R.size();
    int n = B.size();

    std::vector<unsigned int> MAXLENGTH((m)*(m));
    maxlength(R, MAXLENGTH, m, equality);
    vector<char> cigarVector;

    for(int k = 0; k<=m; k++)
    {
        //if (k==0) k=64;
        int ret = algorithm3(R, B, MAXLENGTH, m, n, k, equality, false, cigar, cigarVector);
        if(cigar) {
            cigarOutput = standardCigarOutput(cigarVector);
        }
        if(ret != -1) return ret;
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

string standardCigarOutput(vector<Triple>& Sij)
{
    ostringstream stream;
    int numInsertions = 0, numDeletions = 0, align = 0;
    for(int seqNum = 0; seqNum<Sij.size(); seqNum++) 
    {
        Triple seq = Sij[seqNum];
        if(!seqNum && (seq.f != 0 || seq.c != 0)) { stream << "POS: " << seq.p << endl; align = seq.f; }
        else if(seq.f==0 && seq.c==0)
        {
            if(seqNum) numInsertions++;
            else stream << "POS: " << seq.p << endl;
            align++;
        }
        else
        {
            Triple previous = Sij[seqNum-1-numInsertions < 0 ? 0 : seqNum-1-numInsertions];
            //if (numInsertions!=0 && previous.c+previous.f+numInsertions != seq.c) printf("%dI", numInsertions);
            if(previous.c+previous.f+numInsertions == seq.c && previous.p+previous.f+numInsertions == seq.p)
            {
                numInsertions = 0;
                align += seq.f;
            }
            else if (previous.p+previous.f+numInsertions == seq.p && previous.c+previous.f+numInsertions > seq.c)
            {
                stream << "" << align-numInsertions << "M";
                
                if (numInsertions) stream << "" << numInsertions << "D";
                align = seq.f;
                numInsertions = 0;
            }
            else if (previous.p+previous.f+numInsertions == seq.p && previous.c+previous.f+numInsertions < seq.c)
            {
                if(previous.p+previous.f == Sij[seqNum-1].p) stream << "" << align << "M";
                else stream << "" << align-numInsertions << "M";
                stream << "" << seq.c-(previous.c+previous.f+numInsertions) << "I";
                align = seq.f;
                numInsertions = 0;
            }
            else
            {
                align += seq.f;
                numInsertions = 0;
            }
        }
    }
    if(align) stream << "" << align << "M";
    return stream.str();
}

//TODO
string extendedCigarOutput(vector<Triple>& Sij)
{
    return "";
}

inline int max(int i1, int i2, int i3, int* num)
{
    //printf("%d,%d,%d\n", i1,i2,i3);
    if(i1>=i2 && i1>=i3) {*num = 1; return i1;}
    if(i2>=i1 && i2>=i3) {*num = 2; return i2;}
    *num = 3;
    return i3;
}

bool operator==(const L& lhs ,const L& rhs)
{
        return lhs.d==rhs.d && lhs.e==rhs.e;
}