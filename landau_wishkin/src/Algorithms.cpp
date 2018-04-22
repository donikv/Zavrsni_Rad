#include "Algorithms.hpp"

int algorithm2(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, std::unordered_map<L, int, Hasher, EqualFn>& D, int m, int n, int bStart, int k, int nk, EqualityDefinition& equality)
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
        for(int d = -e; d<=e; d++){
            if(d==-e){
                row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, -5, &num);
            } else if (d==e) {
                row = max(D[L{d,e-1}]+1, -5, D[L{d-1,e-1}], &num);
            } else {
                row = max(D[L{d,e-1}]+1, D[L{d+1,e-1}]+1, D[L{d-1,e-1}], &num);
            }

            while(equality.areEqual(R[row], B[row+d+bStart]) && row<m) row++;
            D[L{d,e}] = row;

            if(row == m){
                printf("\n");
                return e;
            }
        }
    }

    return -1;
}

int algorithm3(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, const std::vector<unsigned int>& MAXLENGTH,int m, int n, int kmin, int kmax, EqualityDefinition& equality, bool prefix, bool cigar, string& cigarOutput)
{
    for(int k = kmin; k<=kmax; k++)
    {
        int j = 0;
        std::vector<Triple> Sij;

        for(int i=0; i<(prefix ? 1 : n-m+k);i++){
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

                        if(f>=1) {
                            if(f != MAXLENGTH[c*m+row]) {
                                row += std::min(f,MAXLENGTH[c*m+row]);
                                goto inst5;
                            } else {
                                row += f;
                            }
                        } else {
                            if(!(equality.areEqual(R[row], B[row+d+i]))){
                                goto inst5;
                            } else {
                                row++;
                            }
                        }
                    }

                    while(equality.areEqual(R[row], B[row+d+i]) && row<m) row++;

                inst5:
                    D[L{d,e}] = row;

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

            for(int l = -k; l<=k; l++){
                std::vector<Triple> sequenceForCurrentL = lSeqMap[L{l,k}];
                if(sequenceForCurrentL.size()<=0) continue;

                if((sequenceForCurrentL.back().p+sequenceForCurrentL.back().f)>=j) {
                    Sij = sequenceForCurrentL;
                    break;
                }
            }
            if(row == m) { if (cigar) cigarOutput = standardCigarOutput(Sij); return k; }
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
    for(int i=4;i>0;i/=2)
    {  
        if(nk != 0) nk = m/nk;
        int k = algorithm2(R,B,D,m,n,0,m/i, nk, equality);
        if (k!=-1) return k;
        nk = i;
    }
    int k = algorithm2(R, B, D, m, n, 0, m/4, 0, equality);
    if (k!=-1) return k;
    return algorithm2(R, B, D, m, n, 0, m, m/2, equality);
}

int findAlgimentWithLowestKGLOBAL(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar, string& cigarOutput)
{
    int m = R.size();
    int n = B.size();

    std::unordered_map<L, int, Hasher, EqualFn> D;

    return (findAlgimentWithLowestKPREFIX(R, B, equality) + n - m);
}

int findAlgimentWithLowestKINFIX(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar, string& cigarOutput)
{
    int m = R.size();
    int n = B.size();

    std::vector<unsigned int> MAXLENGTH((m)*(m));
    maxlength(R, MAXLENGTH, m, equality);

    return algorithm3(R, B, MAXLENGTH, m, n, 0, m, equality, true, cigar, cigarOutput);
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