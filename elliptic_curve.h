#pragma once
#include <iostream>
#include <gmp.h>
#include <vector>
#include <gmpxx.h>


struct _point;

class elliptic_curve {
public:
    struct _point {
        friend class elliptic_curve;
    private:
        mpz_class _x;
        mpz_class _y;
        bool _z = true;

    private:
        _point(const mpz_class&, const mpz_class&, const bool);
    public:
        mpz_class get_x() const;
        mpz_class get_y() const;
        bool get_z() const;
    };

    elliptic_curve(const mpz_class&, const mpz_class&, const mpz_class&);

    _point new_point(const mpz_class&, const mpz_class&, const bool) const;
    _point new_point(const long int, const long int, const bool) const;

    _point generate_point();
    _point generate_point(mpz_class);

    _point sum(const _point&, const _point&) const;
    _point sub(const _point&, const _point&) const;
    _point neg(const _point&) const;

    // Double and add
    _point mul1(const _point&, const mpz_class&) const;

    // Alg 7.2.4: Add and sub
    _point mul2(const _point&, const mpz_class&) const;

    // With ternary expansion
    _point mul3(const _point&, const mpz_class&) const;


    mpz_class get_a() const;
    mpz_class get_b() const;
    mpz_class get_p() const;

    void set_a(const mpz_class&); 
    void set_b(const mpz_class&);

    friend mpz_class order(const elliptic_curve&);


private:
    _point x2(const _point&) const;
    bool exist_point(const mpz_class&, const mpz_class&, const bool) const;
    static std::vector<mpz_class> cross(std::vector<mpz_class>, std::vector <mpz_class>);
    static mpz_class ind(const std::vector<mpz_class>&, const mpz_class&);
    static std::vector<mpz_class> shanks(const elliptic_curve&, std::vector<mpz_class>&, std::vector<mpz_class>&, const elliptic_curve::_point&, const mpz_class&);
    static mpz_class return_t(const elliptic_curve::_point&, const elliptic_curve&, const mpz_class&, const mpz_class&, const mpz_class&, const mpz_class&);

private:

    mpz_class _a;
    mpz_class _b;
    mpz_class _p;

};

std::ostream& operator<<(std::ostream&, const elliptic_curve::_point&);


typedef struct elliptic_curve::_point point;
