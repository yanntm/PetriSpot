#ifndef RATIONAL_H
#define RATIONAL_H

#include <iostream>

template<typename T>
class Rational {
public:
    T numerator;
    T denominator;

    // Constructor with default denominator of 1
    Rational(T num = 0, T den = 1) : numerator(num), denominator(den) {
        simplify();
    }

    // Simplify the fraction using GCD
    void simplify() {
        if (denominator == 0) {
            throw std::invalid_argument("Denominator cannot be zero.");
        }
        T gcd = std::gcd(numerator, denominator);
        numerator /= gcd;
        denominator /= gcd;
        if (denominator < 0) {  // Ensure denominator is positive
            numerator = -numerator;
            denominator = -denominator;
        }
    }

    // Compare with another Rational
    bool operator==(const Rational& other) const {
        return numerator == other.numerator && denominator == other.denominator;
    }

    // Compare with an integer (no temporary Rational needed)
    bool operator==(T integer) const {
        return denominator == 1 && numerator == integer;
    }

    // Symmetry: allow integer == Rational
    friend bool operator==(T integer, const Rational& r) {
        return r == integer;
    }

    friend std::ostream& operator<< (std::ostream &os, const Rational &r)
    {
      os << r.numerator << "/" << r.denominator;
      return os;
    }
};

#endif
