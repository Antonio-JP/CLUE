#ifndef CLUE_RF_H
#define CLUE_RF_H

#include <map>
#include <list>
#include <string>
#include <boost/rational.hpp>

using namespace std;
using namespace boost;

class Monomial {
    private:
        map<int,int> degrees;

        void init(map<int,int>);
    public:
        Monomial(map<int,int> degrees) { // Usual constructor using the map as the indication of degrees
            this->init(degrees);
        }
        Monomial() { // Constructor to build the "empty" monomial, i.e., the monomial 1.
            this->init(map<int,int>());
        }
        Monomial(vector<int> degrees) { // Consider each element as a degree of the index it belongs
            map<int,int> final_input;
            for (int i = 0; i < degrees.size(); i++) {
                if (degrees[i] != 0) {
                    final_input[i] = degrees[i];
                }
            }
            this->init(final_input);
        }

        string to_string(map<int,string> varnames); // Creates a string that represent the monomial given the variable names
        string to_string(vector<string> varnames) { // Creates a string that represent the monomial given the variable names
            map<int,string> final_input;
            for (int i = 0; i < varnames.size(); i++) {
                final_input[i] = varnames[i];
            }

            return this->to_string(final_input);
        }

        Monomial operator*(const Monomial other); // Compute the multiplication of two monomials
        Monomial gcd(const Monomial other); // Compute the GCD of two monomials
        Monomial lcm(const Monomial other); // Compute the LCM of two monomials
        bool operator==(const Monomial other); // Check equality of two monomials
};

#endif