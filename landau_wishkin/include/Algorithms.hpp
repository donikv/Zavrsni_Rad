#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <vector>
#include <unordered_map>
#include <string>
#include <sstream>
#include "Maxlength.hpp"
#include "EqualityDefinition.hpp"

using namespace std;

static string nullString = "NULL";

struct L {
    int d;
    int e;
};

struct Cigar {
    char type;
    unsigned int num;
};

struct Triple {
    unsigned int p;
    unsigned int c;
    unsigned int f;
};

bool operator==(const L& lhs ,const L& rhs);

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

inline int max(int i1, int i2, int i3, int* num);

string standardCigarOutput(vector<Triple>& Sij);
string extendedCigarOutput(vector<Triple>& Sij);

int algorithm3(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, const std::vector<unsigned int>& MAXLENGTH,int m, int n, int kmin, int kmax, EqualityDefinition& equality, bool prefix, bool cigar=false, string& cigarOutput = nullString);

int algorithm2(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, std::unordered_map<L, int, Hasher, EqualFn>& D, int m, int n, int bStart, int k, int nk, EqualityDefinition& equality);

//prefix
int findAlgimentWithLowestKPREFIX(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar=false, string& cigarOutput = nullString);

//global
int findAlgimentWithLowestKGLOBAL(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar=false, string& cigarOutput = nullString);

//global
int findAlgimentWithLowestKINFIX(const std::vector<unsigned char>& R, const std::vector<unsigned char>& B, EqualityDefinition& equality, bool cigar=false, string& cigarOutput = nullString);

#endif