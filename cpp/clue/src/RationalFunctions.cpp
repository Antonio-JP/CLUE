#include "RationalFunctions.hpp"

/* Generic init method for Monomials */
void Monomial::init(map<int,int> degrees) {
    this->degrees = map<int,int>();

    for (pair<const int,int> ppair : degrees) { // We copy the dictionary to the degrees of the monomial
        if (ppair.second < 0) { throw std::invalid_argument("Negative exponent not valid for monomial"); }
        else if (ppair.second > 0) {
            this->degrees[ppair.first] = ppair.second;
        }
    }
}

/* String methods for Monomials */
string Monomial::to_string(map<int,string> varnames) {
    if (this->degrees.size() == 0) {
        return "1";
    }
    string output = "";

    // Here we know there are at least one
    map<int,int>::iterator ppair = this->degrees.begin();
    output += varnames[ppair->first];
    if (ppair->second > 1) {
        output += "^" + std::to_string(ppair->second);
    }

    ppair = next(ppair);
    while(ppair != this->degrees.end()) {
        output += "*" + varnames[ppair->first];
        if (ppair->second > 1) {
            output += "^" + std::to_string(ppair->second);
        }
        ppair = next(ppair);
    }

    return output;
}

/* Operation methods for Monomials*/
Monomial Monomial::operator*(const Monomial other) {
    map<int,int> final_dict;

    for (pair<const int, int> ppair : this->degrees) {
        final_dict[ppair.first] = ppair.second;
    }

    for (pair<const int, int> ppair : other.degrees) {
        int key = ppair.first;
        int degree = ppair.second;

        map<int,int>::const_iterator search = final_dict.find(key); 
        if (search != final_dict.end()) {
            final_dict[key] = search->second + degree;
        } else {
            final_dict[key] = degree;
        }
    }

    return Monomial(final_dict);
}

Monomial Monomial::gcd(const Monomial other) {
    map<int,int> final_dict;

    for (pair<const int, int> this_pair : this->degrees) {
        int key = this_pair.first;
        int degree = this_pair.second;

        map<int,int>::const_iterator search = other.degrees.find(key);
        if (search != other.degrees.end()) {
            final_dict[key] = min(degree, search->second);
        }
    }

    return Monomial(final_dict);
}

Monomial Monomial::lcm(const Monomial other) {
    map<int,int> final_dict;

    for (pair<const int, int> ppair : this->degrees) {
        final_dict[ppair.first] = ppair.second;
    }

    for (pair<const int, int> ppair : other.degrees) {
        int key = ppair.first;
        int degree = ppair.second;

        map<int,int>::const_iterator search = final_dict.find(key);
        if (search != final_dict.end()) {
            final_dict[key] = max(search->second, degree);
        } else {
            final_dict[key] = degree;
        }
    }

    return Monomial(final_dict);
}

bool Monomial::operator==(const Monomial other) {
    for (pair<const int, int> ppair : this->degrees) {
        map<int,int>::const_iterator search = other.degrees.find(ppair.first);
        if (search == other.degrees.end() || search->second != ppair.second) { // Variable not in other or different degree
            return false;
        }
    }
    return true; // all checks were successful
}

int Monomial::degree() {
    int degree = 0;
    for (pair<const int, int> ppair : this->degrees) {
        degree += ppair.second;
    }

    return degree;    
}

int Monomial::degree(int variable) {
    map<int,int>::const_iterator search = this->degrees.find(variable);
    if (search != this->degrees.end()) {
        return search->second;
    }
    return 0;
}