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

        int degree(); // Compute the degree of a monomial
        int degree(int variable); // Compute the degree of a monomial wrt a variable
};

template <typename T>
class SparsePolynomial {
    private:
        map<Monomial,T> coeffs;
        vector<string> varnames;

    public:
        /*********************************************************************/
        /* CONSTRUCTORS */
        SparsePolynomial(vector<string> varnames); // Creates the polynomial 0 with given variables
        SparsePolynomial(vector<string> varnames, map<Monomial,T> data); // Creates the polynomial with existing data
        SparsePolynomial(vector<string> varnames, vector<T> data); // Creates a linear polynomial from  vector of coefficients
        SparsePolynomial(vector<string> varnames, T data); // Creates a constant polynomial
        SparsePolynomial(vector<string> varnames, string data); // Creates the polynomial from a string representation

        /*********************************************************************/
        /* ATTRIBUTE AND PROPERTIES */
        int size(); // Get number of non-zero terms
        vector<string> gens(); // Return the variable names of the polynomial
        vector<Monomial> monomials(); // Return a list with monomials in the polynomial
        vector<T> coefficients(); // Return a list with the coefficients of the polynomial. Same order as :func:`monomials`.
        T content(); // Return the content of the polynomial (GCD of the coefficients)
        T constant_term(); // Return the coefficient of the monomial 1. Return 0 if not in the polynomial.
        T ct() { return this->constant_term(); } // Alias for :func:`constant_term`.
        int degree() { // Return the degree of the polynomial as the maximal degree on the monomials
            int degree = -1;
            for (Monomial mon : this->monomials()) {
                int mon_deg = mon.degree();
                if (mon_deg > degree) { degree = mon_deg; }
            }
            return degree;
        }
        int degree(int variable) { // Return the degree of the polynomial as the maximal degree on the monomials
            int degree = -1;
            for (Monomial mon : this->monomials()) {
                int mon_deg = mon.degree(variable);
                if (mon_deg > degree) { degree = mon_deg; }
            }
            return degree;
        }
        int degree(string variable) { // Return the degree of the polynomial as the maximal degree on the monomials
            int variable_index = -1;
            for (int i = 0; i < this->varnames.size(); i++) {
                if (this->varnames[i] == variable) {
                    variable_index = i; break;
                }
            }
            if (variable_index == -1) {
                throw std::invalid_argument("Variable '" + variable + "' not found in polynomial.");
            }
            return this->degree(variable_index);
        }
        vector<string> variables(); // Names of variables that appear in the polynomial
        
        /*********************************************************************/
        /* ARITHMETIC METHODS */
        SparsePolynomial<T> operator+(SparsePolynomial<T>);
        SparsePolynomial<T> operator+(T);
        SparsePolynomial<T> operator+=(SparsePolynomial<T>);
        SparsePolynomial<T> operator+=(T);
        SparsePolynomial<T> operator-(SparsePolynomial<T>);
        SparsePolynomial<T> operator-(T);
        SparsePolynomial<T> operator-=(SparsePolynomial<T>);
        SparsePolynomial<T> operator-=(T);
        SparsePolynomial<T> operator*(SparsePolynomial<T>);
        SparsePolynomial<T> operator*(T);
        SparsePolynomial<T> operator*=(SparsePolynomial<T>);
        SparsePolynomial<T> operator*=(T);
        SparsePolynomial<T> operator^(int);
        SparsePolynomial<T> operator^=(int);
        SparsePolynomial<T> operator==(SparsePolynomial<T>);
        SparsePolynomial<T> operator==(T);
        SparsePolynomial<T> operator%(SparsePolynomial<T>);
        SparsePolynomial<T> operator%(T);
        SparsePolynomial<T> operator%=(SparsePolynomial<T>);
        SparsePolynomial<T> operator%=(T);
        T operator[](Monomial);
        SparsePolynomial<T> operator()(map<string,SparsePolynomial<T>>); // Substitues variables in the arguments
        SparsePolynomial<T> operator()(vector<SparsePolynomial<T>>); // Substitues variables in the arguments: all variables must be given
        
        /*********************************************************************/
        /* BOOLEAN METHODS */
        bool is_zero() { return (this->size() == 0); } // Checks if the polynomial is the zero polynomial
        bool is_one() { return (this->degree() == 0 && this->ct() == ((T)1)); } // Checks if the polynomial is the 1 polynomial
        bool is_constant() { return (this->degree() == 0); } // Checks if the polynomial is a constant
        bool is_linear() { return (this->degree() == 1 && this->ct() == ((T) 0)); } // Checks if the polynomial is linear
        
        /*********************************************************************/
        /* REPRESENTATION METHODS */
        string to_string(); // Generates a string representation of a polynomial

        /*********************************************************************/
        /* OTHER METHODS */
        vector<T> automated_diff(vector<T>); // Computes the automated differentiation of a polynomial in a given point
        SparsePolynomial<T> derivative(); //Computes the formal derivative of a polynomial
};

template <typename T>
class RationalFunction {
    private:
        SparsePolynomial<T> num, denom;

        void simplify(); // Simplify the numerator and denominator as much as possible
    public:
        /*********************************************************************/
        /* CONSTRUCTORS*/
        RationalFunction(vector<string> varnames, SparsePolynomial<T> num, SparsePolynomial<T> denom);
        RationalFunction(vector<string> varnames, T value): RationalFunction(varnames, SparsePolynomial<T>(varnames, value), SparsePolynomial<T>(varnames, ((T) 1))) { } // Creates the rational function that is constant
        RationalFunction(vector<string> varnames) : RationalFunction(varnames, ((T) 0)) { } // Creates the rational function 0
        RationalFunction(vector<string> varnames, string data);

        /*********************************************************************/
        /* ATTRIBUTES AND PROPERTIES */
        int size() { return this->num.size() + this->denom.size(); }
        vector<string> gens() { return this->num.gens(); }
        T constant_term() { return this->num.ct() / this->denom.ct(); }
        T ct() { return this->constant_term(); } // Alias for :func:`constant_term`.
        int valuation(int variable) { // Return the valuation of a rational function wrt a variable
            return this->num.degree(variable) - this->denom.degree(variable);
        }
        int valuation(string variable) { // Return the valuation of a rational function wrt a variable (with a string input)
            return this->num.degree(variable) - this->denom.degree(variable);
        }
        vector<string> variables() { // Return a list of variables that actually appear in the rational function
            vector<string> result, v1 = this->num.variables(), v2 = this->denom.variables();
            sort(v1.begin(), v1.end());
            sort(v2.begin(), v2.end());
            set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(result));
            return result;
        }
        SparsePolynomial<T> get_poly() {
            if (this->is_polynomial()) {
                return this->num * (((T) 1)/this->denom.ct());
            }
            throw std::domain_error("Polynomial only possible from constant denominators");
        }

        /*********************************************************************/
        /* ARITHMETIC METHODS */
        RationalFunction<T> operator+(RationalFunction<T>);
        RationalFunction<T> operator+(SparsePolynomial<T>);
        RationalFunction<T> operator+(T);
        RationalFunction<T> operator+=(RationalFunction<T>);
        RationalFunction<T> operator+=(SparsePolynomial<T>);
        RationalFunction<T> operator+=(T);
        RationalFunction<T> operator-(RationalFunction<T>);
        RationalFunction<T> operator-(SparsePolynomial<T>);
        RationalFunction<T> operator-(T);
        RationalFunction<T> operator-=(RationalFunction<T>);
        RationalFunction<T> operator-=(SparsePolynomial<T>);
        RationalFunction<T> operator-=(T);
        RationalFunction<T> operator*(RationalFunction<T>);
        RationalFunction<T> operator*(SparsePolynomial<T>);
        RationalFunction<T> operator*(T);
        RationalFunction<T> operator*=(RationalFunction<T>);
        RationalFunction<T> operator*=(SparsePolynomial<T>);
        RationalFunction<T> operator*=(T);
        RationalFunction<T> operator^(int);
        RationalFunction<T> operator^=(int);
        RationalFunction<T> operator==(RationalFunction<T>);
        RationalFunction<T> operator==(SparsePolynomial<T>);
        RationalFunction<T> operator==(T);
        RationalFunction<T> operator/(RationalFunction<T>);
        RationalFunction<T> operator/(SparsePolynomial<T>);
        RationalFunction<T> operator/(T);
        RationalFunction<T> operator/=(RationalFunction<T>);
        RationalFunction<T> operator/=(SparsePolynomial<T>);
        RationalFunction<T> operator/=(T);
        RationalFunction<T> operator()(map<string,SparsePolynomial<T>>); // Substitues variables in the arguments
        RationalFunction<T> operator()(vector<SparsePolynomial<T>>); // Substitues variables in the arguments: all variables must be given

        /*********************************************************************/
        /* BOOLEAN METHODS */
        bool is_polynomial() { return this->denom.is_constant(); }
        bool is_zero() { return this->num.is_zero(); }
        bool is_constant() { return (this->is_zero() || (this->num.is_constant() && this->denom.is_constant())); }

        /*********************************************************************/
        /* REPRESENTATION METHODS */
        string to_string() { // Generates a string representation of a polynomial
            return "(" + this->num.to_string() + ") / (" + this->denom.to_string() + ")";
        }
        
        /*********************************************************************/
        /* OTHER METHODS */
        vector<T> automated_diff(vector<T>); // Computes the automated differentiation of a rational function in a given point
        RationalFunction<T> derivative() { //Computes the formal derivative of a rational function
            SparsePolynomial<T> dnum = this->num.derivative(), ddenom = this->denom.derivative(); 

            return RationalFunction(this->gens(), this->num * ddenom - this->denom * dnum, this->denom^2);
        }
};

#endif