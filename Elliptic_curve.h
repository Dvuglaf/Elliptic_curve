#pragma once
#include <iostream>
#include <gmp.h>
#include <vector>
#include <gmpxx.h>


struct _Point;

class elliptic_curve {
public:
    struct _Point {
        friend class elliptic_curve;
    private:
        mpz_class _x;
        mpz_class _y;
        bool _z = true;

    private:
        _Point(const mpz_class&, const mpz_class&, const bool);
    public:
        mpz_class get_x() const;
        mpz_class get_y() const;
        bool get_z() const;
    };

    elliptic_curve(const mpz_class&, const mpz_class&, const mpz_class&);

    _Point new_point(const mpz_class&, const mpz_class&, const bool) const;
    _Point new_point(const long int, const long int, const bool) const;

    _Point generate_point();
    _Point generate_point(mpz_class);

    _Point sum(const _Point&, const _Point&) const;
    _Point sub(const _Point&, const _Point&) const;
    _Point neg(const _Point&) const;

    // Double and add
    _Point mul1(const _Point&, const mpz_class&) const;

    // Alg 7.2.4: Add and sub
    _Point mul2(const _Point&, const mpz_class&) const;

    // With ternary expansion
    _Point mul3(const _Point&, const mpz_class&) const;


    mpz_class get_a() const;
    mpz_class get_b() const;
    mpz_class get_p() const;

    void set_a(const mpz_class&); 
    void set_b(const mpz_class&);

    friend mpz_class order(const elliptic_curve&);


private:
    _Point x2(const _Point&) const;
    bool exist_point(const mpz_class&, const mpz_class&, const bool) const;
    static std::vector<mpz_class> cross(std::vector<mpz_class>, std::vector <mpz_class>);
    static mpz_class ind(const std::vector<mpz_class>&, const mpz_class&);
    static std::vector<mpz_class> shanks(const elliptic_curve&, std::vector<mpz_class>&, std::vector<mpz_class>&, const elliptic_curve::_Point&, const mpz_class&);
    static mpz_class return_t(const elliptic_curve::_Point&, const elliptic_curve&, const mpz_class&, const mpz_class&, const mpz_class&, const mpz_class&);

private:

    mpz_class _a;
    mpz_class _b;
    mpz_class _p;

};

std::ostream& operator<<(std::ostream&, const elliptic_curve::_Point&);


typedef struct elliptic_curve::_Point point;
