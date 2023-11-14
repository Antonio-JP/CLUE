#include <iostream>
#include <map>
#include <vector>
#include "RationalFunctions.h"

void Monomial_creation() {
    map<int,int> map_in = {{0, 5}, {2, 1}};
    vector<int> vector_in = {5, 4, 3, 0, 1};
    Monomial mon1, mon2, mon3;

    vector<string> varnames = {"x", "y", "z", "a", "b"};

    mon1 = Monomial(); // Must be 1
    mon2 = Monomial(vector_in); // Must be "x^5*y^4*z^3*b"
    mon3 = Monomial(map_in); // Must be "x^5*z"

    if (mon1.to_string(varnames) != "1") { throw std::runtime_error("Error in Monomial creation: the '1' is not achieved: " + mon1.to_string(varnames)); }
    if (mon2.to_string(varnames) != "x^5*y^4*z^3*b") { throw std::runtime_error("Error in Monomial creation: the vector as input did not work: " + mon2.to_string(varnames)); }
    if (mon3.to_string(varnames) != "x^5*z") { throw std::runtime_error("Error in Monomial creation: the dictionary as input did not work: " + mon3.to_string(varnames)); }
}

void Monomial_creation_errors() {
    map<int,int> map_in = {{0, -5}, {2, 1}};
    vector<int> vector_in = {5, 4, -3, 0, 1};
    Monomial mon;

    try {
        mon = Monomial(map_in);
    } catch (std::invalid_argument) {
        try {
            mon = Monomial(vector_in);
        } catch (std::invalid_argument) {
            return;
        }
        throw std::runtime_error("Error in Monomial errors: negative exponent not detected with vector input.");
    }
    throw std::runtime_error("Error in Monomial errors: negative exponent not detected with map input.");
}

void Monomial_operation() {
    // Testing the multiplication, GCD and LCM of monomials
    Monomial mon1 = Monomial({2,5,0,0,1}), mon2 = Monomial({1,7,3,0,0});
    vector<string> varnames = {"x", "y", "z", "a", "b"};
    
    Monomial prod = mon1 * mon2;
    Monomial gcd = mon1.gcd(mon2);
    Monomial lcm = mon1.lcm(mon2);

    if (prod.to_string(varnames) != "x^3*y^12*z^3*b") { throw std::runtime_error("Error in product of Monomials: " + prod.to_string(varnames)); }
    if (gcd.to_string(varnames) != "x*y^5") { throw std::runtime_error("Error in GCD of Monomials: " + gcd.to_string(varnames)); }
    if (lcm.to_string(varnames) != "x^2*y^7*z^3*b") { throw std::runtime_error("Error in LCM of Monomials: " + lcm.to_string(varnames)); }
}

void Monomial_equality() {
    // Testing the multiplication, GCD and LCM of monomials
    Monomial mon1 = Monomial({2,5,0,0,1}), mon2 = Monomial({1,7,3,0,0});
    vector<string> varnames = {"x", "y", "z", "a", "b"};
    
    Monomial prod = mon1 * mon2;
    Monomial gcd = mon1.gcd(mon2);
    Monomial lcm = mon1.lcm(mon2);

    if (!(prod == gcd * lcm)) {throw std::runtime_error("Error in equality: " + prod.to_string(varnames) + " vs " + (gcd * lcm).to_string(varnames)); }
}

void Monomial_degree() {
    // Testing the method degrees
    Monomial mon1 = Monomial({2,5,0,0,1}), mon2 = Monomial({1,7,3,0,0});

    if (mon1.degree() != 8) { throw std::runtime_error("Error in degrees: expected 8 but got " + to_string(mon1.degree())); }
    if (mon2.degree() != 11) { throw std::runtime_error("Error in degrees: expected 11 but got " + to_string(mon2.degree())); }
    if (mon1.degree(0) != 2) { throw std::runtime_error("Error in degrees with variable: expected 2 but got " + to_string(mon1.degree(0))); }
    if (mon1.degree(1) != 5) { throw std::runtime_error("Error in degrees with variable: expected 5 but got " + to_string(mon1.degree(1))); }
    if (mon1.degree(2) != 0) { throw std::runtime_error("Error in degrees with variable: expected 0 but got " + to_string(mon1.degree(2))); }
    if (mon1.degree(10) != 0) { throw std::runtime_error("Error in degrees with variable: expected 0 but got " + to_string(mon1.degree(10))); }
}

int main() {
    Monomial_creation();
    Monomial_creation_errors();
    Monomial_operation();
    Monomial_equality();
    Monomial_degree();

    return 0;
}