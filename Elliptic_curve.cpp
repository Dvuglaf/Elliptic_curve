// Elliptic_curve.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Elliptic_curve.h"
#include "utils.h"
#include <algorithm>
#include <gmpxx.h>

const mpz_t* point::get_x() const {
    return &_x;
}

const mpz_t* point::get_y() const {
    return &_y;
}

const bool point::get_z() const {
    return _z;
}

point::point(const mpz_t x, const mpz_t y, elliptic_curve& curve, const bool z = true) {
    mpz_inits(_x, _y, NULL);
    mpz_set(_x, x);
    mpz_set(_y, y);
    _curve = &curve;
    _z = z;
}

point::point(const int x, const int y, elliptic_curve& curve, const bool z = true) {
    mpz_inits(_x, _y, NULL);
    mpz_set_si(_x, x);
    mpz_set_si(_y, y);
    _curve = &curve;
    _z = z;
}

elliptic_curve::elliptic_curve(const mpz_t a, const mpz_t b, const mpz_t p) {
    mpz_inits(_a, _b, _p, NULL);
    mpz_set(_a, a);
    mpz_set(_b, b);
    mpz_set(_p, p);
}

void elliptic_curve::set_a(const mpz_t a) {
    mpz_set(_a, a);
}
void elliptic_curve::set_b(const mpz_t b) {
    mpz_set(_b, b);
}
void elliptic_curve::set_p(const mpz_t p) {
    mpz_set(_p, p);
}

bool elliptic_curve::exist_point(const mpz_t x, const mpz_t y, const bool z = true) {
    if (!z)
        return true;
    mpz_t copy_x, copy_y;
    mpz_inits(copy_x, copy_y, NULL);
    mpz_set(copy_x, x);
    mpz_set(copy_y, y);

    mpz_powm_ui(copy_y, copy_y, 2, _p); // y^2
    mpz_t temp;
    mpz_init(temp);
    mpz_powm_ui(temp, copy_x, 3, _p); 
    mpz_mul(copy_x, copy_x, _a); 
    mpz_add(temp, temp, copy_x); 
    mpz_add(temp, temp, _b); 
    mpz_mod(temp, temp, _p); // temp = x^3 + x*a + b mod(_p)
    bool result = static_cast<bool>(mpz_cmp(temp, copy_y));
    mpz_clears(temp, copy_x, copy_y, NULL);
    if (result)
        std::cout << "Point " << point(x, y, *this) << "are not on the curve.\n";
    return !result;
}

point elliptic_curve::new_point(const mpz_t x, const mpz_t y, const bool z = true) {
    if (z == false) {
        return point(x, y, *this, false);
    }
    mpz_t copy_x, copy_y;
    mpz_inits(copy_x, copy_y, NULL);

    mpz_set(copy_x, x);
    mpz_set(copy_y, y);

    mpz_mod(copy_x, copy_x, _p);
    mpz_mod(copy_y, copy_y, _p);

    if (elliptic_curve::exist_point(copy_x, copy_y)) {
        point result(copy_x, copy_y, *this);
        mpz_clears(copy_x, copy_y, NULL);
        return result;
    }
    mpz_clears(copy_x, copy_y, NULL);
    return point(0, 1, *this, false);
}

point elliptic_curve::new_point(long int x, long int y, const bool z = true) {
    if (z == false) {
        return point(x, y, *this, false);
    }
    mpz_t big_x, big_y;
    mpz_inits(big_x, big_y, NULL);
    mpz_set_si(big_x, x);
    mpz_set_si(big_y, y);

    point result = new_point(big_x, big_y);
    mpz_clears(big_x, big_y, NULL);
    return result;
}

point elliptic_curve::sum(const point& P, const point& Q) {
    if (P._z == 0) return Q; // P = O
    if (Q._z == 0) return P; // Q = O

    mpz_t m;
    mpz_init(m);

    if (mpz_cmp(P._x, Q._x) == 0) {

        mpz_t temp;
        mpz_init(temp);
        mpz_add(temp, P._y, Q._y);
        mpz_mod(temp, temp, _p); // temp = y1 + y2 mod(_p)
        if (mpz_cmp_ui(temp, 0) == 0) {
            mpz_clears(m, temp, NULL);
            return point(0, 1, *this, 0);
        }

        mpz_t first;
        mpz_t second;
        mpz_inits(first, second, NULL);

        mpz_powm_ui(first, P._x, 2, _p);
        mpz_mul_ui(first, first, 3);
        mpz_add(first, first, _a); // first = 3(x1)^2 + a

        mpz_mul_si(second, P._y, 2);
        mpz_invert(second, second, _p); // second = (2y1)-^-1 mod(_p)

        mpz_mul(m, first, second);
        mpz_mod(m, m, _p); // m = first * second
        mpz_mod(m, m, _p);

        mpz_clears(first, second, temp, NULL);
    }
    else {
        mpz_t first;
        mpz_t second;
        mpz_inits(first, second, NULL);

        mpz_sub(first, Q._y, P._y); // first = y2 - y1

        mpz_sub(second, Q._x, P._x);
        mpz_invert(second, second, _p); // second = (x2 - x1)^-1

        mpz_mul(m, first, second); // m = first * second
        mpz_mod(m, m, _p);
        mpz_clears(first, second, NULL);
    }

    mpz_t x3;
    mpz_t y3;
    mpz_t temp;
    mpz_inits(x3, y3, temp, NULL);

    mpz_powm_ui(temp, m, 2, _p);
    mpz_sub(x3, temp, P._x); // x3 = m^2 - x1
    mpz_sub(x3, x3, Q._x); 
    mpz_mod(x3, x3, _p); // x3 = m^2 - x1 -x2 mod(_p) 

    mpz_sub(temp, P._x, x3); // temp = x1 - x3
    mpz_mul(y3, temp, m); // y3 = m * temp
    mpz_sub(y3, y3, P._y);
    mpz_mod(y3, y3, _p); // y3 = m * temp - y1 mod(_p)

    point res(x3, y3, *this, 1);
    mpz_clears(x3, y3, temp, m, NULL);

    return res;
}

std::ostream& operator<<(std::ostream& out, const point& P) {
    if (P.get_z() == false) {
        out << "O\n";
    }
    else {
        out << "( " << *P.get_x() << ", " << *P.get_y() << " )";
    }
    return out;
}

point elliptic_curve::x2(const point& P) {
    return sum(P, P);
}

point elliptic_curve::mul(const point& P, const mpz_t n) {
    if (P._z == false) {
        return P;
    }
    mpz_t abs_n;
    mpz_init(abs_n);
    mpz_abs(abs_n, n);
    mpz_t zero;
    mpz_init_set_si(zero, 0);
    bool negative = false;
    if (mpz_cmp(n, zero) == 0) { // 0 * P = 0
        mpz_clears(zero, abs_n, NULL);
        return point(0, 1, *this, false);
    }
    else if (mpz_cmp(n, zero) < 0) { // -n*P = (|n|P)^-1        
        negative = true;
    }
    point result(0, 1, *this, false);
    auto temp = new_point(P._x, P._y);
    auto bits = utils::binary(abs_n); // inverse bits of n
    for (auto& bit : bits) {
        if (bit == '1') {
            result = sum(result, temp); // result = result + temp
        }
        temp = x2(temp); // temp = temp + temp
    }
    if (negative) {
        result = neg(result); // (|n|P)^-1
    }
    mpz_clear(abs_n);
    mpz_mod(result._y, result._y, _p);
    return result;
}

point elliptic_curve::mul(const point& P, const long int n) {
    mpz_t temp;
    mpz_init(temp);
    mpz_set_si(temp, n);
    auto result = mul(P, temp);
    mpz_clear(temp);
    return result;
}

point elliptic_curve::neg(const point& P) {

    mpz_t zero;
    mpz_init_set_si(zero, 0);

    if (P._z == false) {
        mpz_clear(zero);
        return P;
    }
    else if (mpz_cmp(P._y, zero) != 0) {
        mpz_t neg_y;
        mpz_init(neg_y);
        mpz_neg(neg_y, P._y); // neg = -P.y
        point result(P._x, neg_y, *this);
        mpz_clears(neg_y, zero, NULL);
        mpz_mod(result._y, result._y, _p);
        return result;
    }
    else {
        mpz_clear(zero);
        return P;
    }
}

point elliptic_curve::generate_point() {

    gmp_randstate_t rand_state;
    gmp_randinit_mt(rand_state);

    mpz_t random_x;
    mpz_t t;
    mpz_inits(random_x, t, NULL);

    gmp_randseed_ui(rand_state, utils::get_seed());
    
    while (true) {
        mpz_urandomm(random_x, rand_state, _p);
        mpz_powm_ui(t, random_x, 2, _p);
        mpz_add(t, t, _a);
        mpz_mul(t, t, random_x);
        mpz_add(t, t, _b);
        mpz_mod(t, t, _p); //t = (x(x^2 + a) + b (mod p)

        if (mpz_legendre(t, _p) == -1) {
            continue;
        }
        else
            break;
    }
    utils::sqrtm(t, t, _p); // t = sqrt(t) in mod (_p)
    point result(random_x, t, *this);
    mpz_clears(random_x, t, NULL);
    gmp_randclear(rand_state);
    return result;
}

point elliptic_curve::generate_point(const mpz_t x, const elliptic_curve& E) {
    
    mpz_t t;
    mpz_init(t);

    mpz_powm_ui(t, x, 2, E._p);
    mpz_add(t, t, E._a);
    mpz_mul(t, t, x);
    mpz_add(t, t, E._b);
    mpz_mod(t, t, E._p); //t = (x(x^2 + a) + b (mod p)

    if (mpz_legendre(t, E._p) == -1) {
        mpz_clear(t);
        return point (0, 1, *this, false);
    }

    utils::sqrtm(t, t, E._p); // t = sqrt(t) in mod (_p)
    point result(x, t, *this);
    mpz_clear(t);
    return result;
}

std::deque<mpz_class> cross(std::deque<mpz_class> A, std::deque <mpz_class> B) {
    std::cout << "A = {";
    for (auto& it : A) {
        std::cout << it << ", ";
    }
    std::cout << "}\n";
    std::cout << "B = {";
    for (auto& it : B) {
        std::cout << it << ", ";
    }
    std::cout << "}\n";
    std::sort(A.begin(), A.end(), [](mpz_class first, mpz_class second) {
        if (mpz_cmp(first.get_mpz_t(), second.get_mpz_t()) < 0)
            return true;
        return false;
        });

    std::sort(B.begin(), B.end());


    unsigned long i = 0, j = 0;
    std::deque<mpz_class> S;

    std::cout << "Starting cycle in cross...\n";
    while ((i < A.size()) && (j < B.size())) {
        if (A[i] <= B[j]) {
            if (A[i] == B[j]) {
                std::cout << A[i] << " == " << B[j] << std::endl;
                S.push_back(A[i]);
            }
            ++i;
            while ((i < A.size() - 1) && (A[i] == A[i - 1])) {
                ++i;
            }
        }
        else {
            ++j;
            while ((j < B.size() - 1) && (B[j] == B[j - 1])) {
                ++j;
            }
        }
    }
    std::cout << "End cycle in cross\n";
    return S;

}

std::deque<mpz_class> elliptic_curve::shanks(
    std::deque<mpz_class>& A, std::deque<mpz_class>& B, 
    const point P, elliptic_curve E, const mpz_class W, const mpz_class x
) {
    mpz_t idx, bound, temp;
    mpz_inits(bound, temp, NULL);
    mpz_sub_ui(bound, W.get_mpz_t(), 1);
    mpz_init_set_ui(idx, 0);
    for (; mpz_cmp(idx, bound) <= 0; mpz_add_ui(idx, idx, 1)) {
        mpz_set(temp, E._p);
        mpz_add(temp, temp, idx);
        mpz_add_ui(temp, temp, 1);
        point point_temp(0, 1, E, false);
        point_temp = E.mul(P, temp);
        A.push_back(mpz_class(point_temp._x));
    }

    mpz_set_ui(idx, 0);
    mpz_set(bound, W.get_mpz_t());
    for (; mpz_cmp(idx, bound) <= 0; mpz_add_ui(idx, idx, 1)) {
        mpz_set(temp, W.get_mpz_t());
        mpz_mul(temp, temp, idx);
        point point_temp(0, 1, E, false);
        point_temp = E.mul(P, temp);
        B.push_back(mpz_class(point_temp._x));
    }

    mpz_clears(idx, bound, temp, NULL);
    return cross(A, B);
}

mpz_class ind(std::deque<mpz_class>& S, mpz_class s) {
    mpz_t i;
    mpz_init_set_ui(i, 0);
    for (auto& it : S) {
        if (it == s) {
            mpz_clear(i);
            return mpz_class(i);
        }
        mpz_add_ui(i, i, 1);
    }
    mpz_clear(i);
    return mpz_class(-1);

}

mpz_class return_t(const point& P, elliptic_curve E, const mpz_class& beta, const mpz_class& gamma, const mpz_class& W, const mpz_class p) {
    mpz_class temp = beta;
    temp = temp + gamma * W;
    mpz_class test = p + 1 + temp;
    auto pnt = E.mul(P, test.get_mpz_t());
    if (pnt.get_z() == false)
        return temp;
    else {
        return beta - gamma * W;
    }
}

mpz_class elliptic_curve::order() {
    if (mpz_cmp_ui(_p, 229) <= 0) {

    }

    elliptic_curve E(_a, _b, _p);

    mpz_t g;
    mpz_init(g);
    utils::rand_not_sqr_res(g, _p);

    mpz_t W;
    mpz_init(W);
    mpf_t temp, sqrt_2;
    mpf_inits(temp, sqrt_2, NULL);
    mpf_set_z(temp, _p);
    mpf_sqrt(temp, temp);
    mpf_sqrt(temp, temp);
    mpf_sqrt_ui(sqrt_2, 2);
    mpf_mul(temp, temp, sqrt_2);
    mpf_ceil(temp, temp);
    mpz_set_f(W, temp);
    std::cout << "W = " << W << std::endl;
    mpf_clears(temp, sqrt_2, NULL);

    mpz_t c, d;
    mpz_inits(c, d, NULL);
    mpz_powm_ui(c, g, 2, _p);
    mpz_mul(c, c, _a);
    mpz_mod(c, c, _p);
    
    mpz_powm_ui(d, g, 3, _p);
    mpz_mul(d, d, _b);
    mpz_mod(d, d, _p);
    
    mpz_t sigma, y_2, a_x;
    mpz_inits(sigma, y_2, a_x, NULL);
    gmp_randstate_t rand_state;
    gmp_randinit_mt(rand_state);
    while (true) {

        mpz_t x;
        mpz_init(x);
        gmp_randseed_ui(rand_state, utils::get_seed());
        mpz_urandomm(x, rand_state, _p);


        mpz_powm_ui(y_2, x, 3, _p);
        mpz_set(a_x, _a);
        mpz_mul(a_x, a_x, x);

        mpz_add(y_2, y_2, a_x);
        mpz_add(y_2, y_2, _b);
        mpz_mod(y_2, y_2, _p);

        mpz_set_si(sigma, mpz_legendre(y_2, _p));
        std::cout << "sigma = " << sigma << std::endl;

        if (mpz_cmp_ui(sigma, 0) == 0) {
            std::cout << "legandre(y^2, p) = 0, continue!\n" << std::endl;
            mpz_clear(x);
            continue;
        }
        if (mpz_cmp_ui(sigma, 1) == 0) {
            ;
        }
        else {
            std::cout << "set (c, d) = (" << c << ", " << d << " )" <<std::endl;
            E.set_a(c);
            E.set_b(d);
            mpz_mul(x, x, g);
            std::cout << "g = " << g << std::endl;
            std::cout << "set x = " << x << std::endl;
        }
        std::cout << "Starting generate point..." << std::endl;
        point P = generate_point(x, E);
        std::cout << "Success.\n";
        if (P._z == false) {
            std::cout << "for this x point does not exist, continue!\n\n";
            mpz_clear(x);
            E.set_a(_a);
            E.set_b(_b);
            continue;
        }
        std::cout << "point P : " << P << std::endl;

        std::deque<mpz_class> A, B, S;

        S = shanks(A, B, P, E, mpz_class(W), mpz_class(x));
        std::cout << "S = {";
        for (auto& it : S) {
            std::cout << it << ", ";
        }
        std::cout << "}\n";

        if (S.size() != 1) {
            mpz_clear(x);
            E.set_a(_a);
            E.set_b(_b);
            continue;
        }

        gmp_randclear(rand_state);
        auto s = S[0];
        auto beta = ind(A, s);
        auto gamma = ind(B, s);
        std::cout << "beta = " << beta << std::endl;
        std::cout << "gamma = " << gamma << std::endl;

        auto t = return_t(P, E, beta, gamma, mpz_class(W), mpz_class(_p));
        std::cout << "t = " << t << std::endl;
        mpz_t result;
        mpz_init_set(result, sigma);
        mpz_mul(result, t.get_mpz_t(), result);
        mpz_add(result, result, _p);
        mpz_add_ui(result, result, 1);
        mpz_class ret(result);

        mpz_clears(W, y_2, a_x, x, result, NULL);
        return mpz_class(result);
    }

}

int main() {

    mpz_class a = 2, b = 3, p("599234844171323798014378139219684288087362388635279936453684278910360928027555754598727507571581246001777277670122448689713029876360082601445563780429103017907332764522609770882588110158952861810496271751306659999739053379195854324818229379783745932907778328817332573106903755814283051471753305325707");
    mpz_class c = 120, d = 1, e = 327, f = 12421;

    std::deque<mpz_class> vec;
    vec.push_back(b);
    vec.push_back(a);
    vec.push_back(p);
    vec.push_back(c);
    vec.push_back(e);
    vec.push_back(f);
    vec.push_back(d);
    vec.push_back(e);
    std::sort(vec.begin(), vec.end());
    
    elliptic_curve curve(a.get_mpz_t(), b.get_mpz_t(), p.get_mpz_t());
    point first = curve.generate_point();
    point second = curve.generate_point();

    std::cout << curve.sum(first, second);
    std::cout << curve.mul(first, f.get_mpz_t());
    std::cout << curve.neg(first);

    //std::cout << curve.order();
/*
    point pnt1 = curve.new_point(3, 6);
    point pnt2 = curve.new_point(80, -10);
    std::cout << curve.sum(pnt1, pnt2);
    std::cout << curve.mul(pnt1, -2);
    std::cout << curve.neg(pnt1);
*/
    //point pnt = curve.generate_point();
    //std::cout << pnt << std::endl;


}



//mpz_set_str(p, "599234844171323798014378139219684288087362388635279936453684278910360928027555754598727507571581246001777277670122448689713029876360082601445563780429103017907332764522609770882588110158952861810496271751306659999739053379195854324818229379783745932907778328817332573106903755814283051471753305325707", 10);

