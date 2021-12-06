#pragma once
#include <iostream>
#include <gmp.h>

struct point;

class elliptic_curve {
private:
    mpz_t _a;
    mpz_t _b;
    mpz_t _p;

public:
    elliptic_curve(const mpz_t, const mpz_t, const mpz_t);

    point new_point(const mpz_t, const mpz_t, const bool);
    point new_point(long int, long int, const bool);
        
    point sum(const point&, const point&);
    point mul(const point&, const mpz_t);
    point mul(const point&, const long int);
    point neg(const point&);

    point generate_point();

private:
    point x2(const point&);
    bool exist_point(const mpz_t, const mpz_t, const bool);
};

struct point {
    friend class elliptic_curve;
private:
    mpz_t _x;
    mpz_t _y;
    elliptic_curve* _curve = nullptr;
    bool  _z = true;

private:
    point(const mpz_t, const mpz_t, elliptic_curve&, const bool);
    point(const int, const int, elliptic_curve&, const bool );

public:
    const mpz_t* get_x() const;
    const mpz_t* get_y() const;
    const bool get_z() const;
};

std::ostream& operator<<(std::ostream&, const point&);



/*
    1. Сложение +
    2. Обратная точка +
    3. Умножение (удвоение) +
    4. Умножение (тернарное разложение)
    5. Рандомная точка на кривой +
    6. Подсчет точек на кривой

    Сделать, чтоб извне нельзя было что то делать с точками и проверки +

*/
