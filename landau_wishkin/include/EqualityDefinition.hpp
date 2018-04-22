#ifndef EQUALITY_DEFINITION_HPP 
#define EQUALITY_DEFINITION_HPP

#include <vector>
#include <cstring>
#include <string>

using namespace std;

typedef uint64_t Word;
static const int WORD_SIZE = sizeof(Word) * 8; // Size of Word in bits
static const Word WORD_1 = (Word)1;
static const Word HIGH_BIT_MASK = WORD_1 << (WORD_SIZE - 1);  // 100..00
static const int MAX_UCHAR = 255;

typedef struct {
    char first;
    char second;
} EdlibEqualityPair;

class EqualityDefinition {
private:
    bool matrix[MAX_UCHAR + 1][MAX_UCHAR + 1];
public:
    EqualityDefinition(const string& alphabet,
                       const EdlibEqualityPair* additionalEqualities = NULL,
                       const int additionalEqualitiesLength = 0);
    /**
     * @param a  Element from transformed sequence.
     * @param b  Element from transformed sequence.
     * @return True if a and b are defined as equal, false otherwise.
     */
    bool areEqual(unsigned char a, unsigned char b) const;
};

string transformSequences(const char* const queryOriginal, const int queryLength,
                                 const char* const targetOriginal, const int targetLength, vector<unsigned char>& queryTransformed, vector<unsigned char>& targetTransformed);

#endif