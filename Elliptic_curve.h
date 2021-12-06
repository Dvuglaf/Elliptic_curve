#pragma once
#include <iostream>
#include <gmp.h>
#include <deque>
#include <gmpxx.h>

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

    mpz_class order();

    point generate_point();
    point generate_point(const mpz_t, const elliptic_curve&);


private:
    point x2(const point&);
    bool exist_point(const mpz_t, const mpz_t, const bool);
    void set_a(const mpz_t);
    void set_b(const mpz_t);
    void set_p(const mpz_t);
    static std::deque<mpz_class> shanks(std::deque<mpz_class>& A, std::deque<mpz_class>& B, const point P, const elliptic_curve E, const mpz_class W, const mpz_class x);

    mpz_t* get_p();

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
std::deque<mpz_class> cross(std::deque<mpz_class>, std::deque < mpz_class>);



/*
    1. Сложение +
    2. Обратная точка +
    3. Умножение (удвоение) +
    4. Умножение (тернарное разложение)
    5. Рандомная точка на кривой +
    6. Подсчет точек на кривой

    Сделать, чтоб извне нельзя было что то делать с точками и проверки +

*/
