#pragma once
#include <iostream>
#include <gmp.h>
#include <vector>
#include <gmpxx.h>
class _Point;


class elliptic_curve {
public:
    struct _Point {
        friend class elliptic_curve;
    private:
        mpz_class _x;
        mpz_class _y;
        bool  _z = true;

    private:
        _Point(const mpz_class, const mpz_class, const bool);
        _Point(const int, const int, const bool);

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
    _Point mul(const _Point&, const mpz_class&) const;
    _Point mul(const _Point&, const long int) const;
    _Point neg(const _Point&) const;


    mpz_class get_a() const;
    mpz_class get_b() const;
    mpz_class get_p() const;

    void set_a(const mpz_class&); 
    void set_b(const mpz_class&);



private:
    _Point x2(const _Point&) const;
    bool exist_point(const mpz_class&, const mpz_class&, const bool) const;

private:

    mpz_class _a;
    mpz_class _b;
    mpz_class _p;

};


std::vector<mpz_class> shanks(const elliptic_curve&, std::vector<mpz_class>&, std::vector<mpz_class>&, const elliptic_curve::_Point&, const mpz_class&, const mpz_class&);
mpz_class order(const elliptic_curve&);

std::ostream& operator<<(std::ostream&, const elliptic_curve::_Point&);
std::vector<mpz_class> cross(std::vector<mpz_class>, std::vector < mpz_class>);

typedef class elliptic_curve::_Point point;

/*
    1. Сложение +
    2. Обратная точка +
    3. Умножение (удвоение) +
    4. Умножение (тернарное разложение)
    5. Рандомная точка на кривой +
    6. Подсчет точек на кривой

    Сделать, чтоб извне нельзя было что то делать с точками и проверки +

*/
