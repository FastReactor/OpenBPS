#ifndef UNCERTAINTY_H
#define UNCERTAINTY_H
#include <exception>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

namespace openbps {
//==============================================================================
// Uncertainty class description
//==============================================================================
template<typename T>
class Uncertainty
{
private:
    T real; //!> real part of number
    T dev;  //!> deviative part of number
public:
    //----------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    Uncertainty(T n = 0, T s = 0) : real(n), dev(s) {}
    ~Uncertainty() = default;
    //----------------------------------------------------------------------------
    //! Methods and overloading operators
    // u1 + u2
    friend const Uncertainty operator+(
        const Uncertainty& u1, const Uncertainty& u2)
    {
        return Uncertainty(
                u1.real + u2.real,
                u1.dev + u2.dev);
    }
    // u1 - u2
    friend const Uncertainty operator-(
        const Uncertainty& u1, const Uncertainty& u2)
    {
        return Uncertainty(
                u1.real - u2.real,
                u1.dev + u2.dev);
    }
    // u1 * u2
    friend const Uncertainty operator*(
        const Uncertainty& u1, const Uncertainty& u2)
    {
        return Uncertainty(
                u1.real * u2.real,
                u1.dev * u2.real + u2.dev * u1.real);
    }
    // u1 / u2
    friend const Uncertainty operator/(
        const Uncertainty& u1, const Uncertainty& u2)
    {
        if (u2.real == 0) {
            throw 0;
        } else {
            return Uncertainty(
                u1.real / u2.real,
                u1.dev * u2.real + u2.dev / u2.real / u2.real * u1.real);
        }
    }
    // cout << u1
    friend std::ostream& operator<< (std::ostream& out,
                                     const Uncertainty& u1) {
        out << u1.real << " +/- " << u1.dev;}
    // >> u1
    friend std::istream& operator>> (std::istream& is,
                                      Uncertainty& u1) {
        is >> u1.real;
        return is;
    }
    // u1 * T
    const Uncertainty operator*(const T& d) {
        return Uncertainty(real * d, dev * d);}
    // u1 / T
    const Uncertainty operator/(const T& d) {
        return Uncertainty(real / d, dev / d);}

    //! Get the real part
    T Real() const  {return real;}
    //! Get the deviation part
    T Dev () const {return dev;}
    //! Add diviation to a number
    void Adddeviation(Uncertainty& val) {
        this->dev = val.Dev();
    }
    //! Add diviation to a number
    void Adddeviation(const T& val) {
        this->dev = val;
    }

};

// Type definitions of uncertanty with ordinary and double precesion
typedef Uncertainty<double>  udouble;
typedef Uncertainty<float>   ufloat;

//==============================================================================
// Non class methods
//==============================================================================

//! Split one Uncertainty vector on two T vectors
//!
//! \param[in] v uncertainty vector type T
//! \return result of calculation
template <typename T>
std::pair<std::vector<T>, std::vector<T>>
usplit(std::vector<Uncertainty<T>>& v) {

    std::vector<T> u1;
    std::vector<T> u2;
    std::for_each(v.begin(), v.end(), [&u1, &u2](Uncertainty<T>& n) {
        u1.push_back(n.Real());
        u2.push_back(n.Dev());}
    );
    std::pair<std::vector<T>, std::vector<T>> result {std::make_pair(u1, u2)};
    return result;
}

//! Unite two T vectors into one Uncertainty
//!
//! \param[in] u1 Real part of new vector
//! \param[in] u2 Deviation part of new vector
//! \return result Uncertainty vector type T with real u1 and Dev u2 parts
template <typename T>
std::vector<Uncertainty<T>> ujoin(const std::vector<T>& u1,
                                  const std::vector<T>& u2) {

    std::vector<Uncertainty<T>> result;
    if (u1.size() != u2.size()) {
        std::cerr << "Cannot unite two vectors into vec Uncertainty "
                  <<u1.size() << " Differse from " << u2.size() << std::endl;
    } else {
        for (size_t i; i < u1.size(); i++) {
            result.push_back(Uncertainty<T>(u1[i], u2[i]));
        }
    }
    return result;
}

} //namespace openbps
#endif // UNCERTAINTY_H
