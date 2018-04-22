#include "EqualityDefinition.hpp"

using namespace std;

EqualityDefinition::EqualityDefinition(const string& alphabet,
                    const EdlibEqualityPair* additionalEqualities,
                    const int additionalEqualitiesLength) {
    for (int i = 0; i < (int) alphabet.size(); i++) {
        for (int j = 0; j < (int) alphabet.size(); j++) {
            matrix[i][j] = (i == j);
        }
    }
    if (additionalEqualities != NULL) {
        for (int i = 0; i < additionalEqualitiesLength; i++) {
            size_t firstTransformed = alphabet.find(additionalEqualities[i].first);
            size_t secondTransformed = alphabet.find(additionalEqualities[i].second);
            if (firstTransformed != string::npos && secondTransformed != string::npos) {
                matrix[firstTransformed][secondTransformed] = matrix[secondTransformed][firstTransformed] = true;
            }
        }
    }
}

/**
 * @param a  Element from transformed sequence.
 * @param b  Element from transformed sequence.
 * @return True if a and b are defined as equal, false otherwise.
 */
bool EqualityDefinition::areEqual(unsigned char a, unsigned char b) const {
    return matrix[a][b];
}

string transformSequences(const char* const queryOriginal, const int queryLength,
                                const char* const targetOriginal, const int targetLength, vector<unsigned char>& queryTransformed, vector<unsigned char>& targetTransformed) {
    string alphabet = "";

    // Alphabet information, it is constructed on fly while transforming sequences.
    // letterIdx[c] is index of letter c in alphabet.
    unsigned char letterIdx[MAX_UCHAR + 1];
    bool inAlphabet[MAX_UCHAR + 1]; // inAlphabet[c] is true if c is in alphabet
    for (int i = 0; i < MAX_UCHAR + 1; i++) inAlphabet[i] = false;

    for (int i = 0; i < queryLength; i++) {
        unsigned char c = static_cast<unsigned char>(queryOriginal[i]);
        if (!inAlphabet[c]) {
            inAlphabet[c] = true;
            letterIdx[c] = (unsigned char) alphabet.size();
            alphabet += queryOriginal[i];
        }
        queryTransformed[i] = letterIdx[c];
    }
    for (int i = 0; i < targetLength; i++) {
        unsigned char c = static_cast<unsigned char>(targetOriginal[i]);
        if (!inAlphabet[c]) {
            inAlphabet[c] = true;
            letterIdx[c] = (unsigned char) alphabet.size();
            alphabet += targetOriginal[i];
        }
        targetTransformed[i] = letterIdx[c];
    }

    return alphabet;
}
