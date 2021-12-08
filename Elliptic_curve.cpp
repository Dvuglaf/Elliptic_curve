// Elliptic_curve.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Elliptic_curve.h"
#include "utils.h"
#include <algorithm>
#include <gmpxx.h>

const mpz_class elliptic_curve::_Point::get_x() const {
    return _x;
}

mpz_class elliptic_curve::_Point::_Point::get_y() const {
    return _y;
}

const bool elliptic_curve::_Point::_Point::get_z() const {
    return _z;
}

 
elliptic_curve::_Point::_Point(const mpz_class x, const mpz_class y, const bool z = true) : _x(x), _y(y), _z(z) {}

elliptic_curve::_Point::_Point(const int x, const int y, const bool z = true) : _x(x), _y(y), _z(z) {}

elliptic_curve::elliptic_curve(const mpz_class a, const mpz_class b, const mpz_class p) : _a(a), _b(b), _p(p) {}

void elliptic_curve::set_a(const mpz_class a) {
    _a = a;
}
void elliptic_curve::set_b(const mpz_class b) {
    _b = b;
}
void elliptic_curve::set_p(const mpz_class p) {
    _p = p;
}

mpz_class elliptic_curve::get_a() const {
    return _a;
}
mpz_class elliptic_curve::get_b() const {
    return _b;
}
mpz_class elliptic_curve::get_p() const {
    return _p;
}

bool elliptic_curve::exist_point(const mpz_class x, const mpz_class y, const bool z = true) const {
    if (!z)
        return true;

    mpz_class left = y * y;
    mpz_mod(left.get_mpz_t(), left.get_mpz_t(), _p.get_mpz_t());
    mpz_class right = (x * x * x + x * _a + _b) % _p;
    mpz_mod(right.get_mpz_t(), right.get_mpz_t(), _p.get_mpz_t());
    if (left != right) {
        std::cout << "Point = (" << x << ", " << y << ") are not on the curve.\n";
        return false;
    }
    return true;
}

elliptic_curve::_Point elliptic_curve::new_point(const mpz_class x, const mpz_class y, const bool z = true) const {
    if (z == false) {
        return point(0, 1, false);
    }
    if (elliptic_curve::exist_point(x, y, z)) {
        return point(x, y, z);
    }
    return point(0, 1, false);
}

elliptic_curve::_Point elliptic_curve::new_point(long int x, long int y, const bool z = true) const {
    if (z == false) {
        return point(0, 1, false);
    }
    return new_point(mpz_class(x), mpz_class(y), z);
}

elliptic_curve::_Point elliptic_curve::sum(const elliptic_curve::_Point& P, const elliptic_curve::_Point& Q) const {
    auto _x1 = P.get_x();
    auto _y1 = P.get_y();
    auto _x2 = Q.get_x();
    auto _y2 = Q.get_y();

    if (P.get_z() == 0) return Q; // P = O
    if (Q.get_z() == 0) return P; // Q = O

    mpz_class m;

    if (_x1 == _x2) {
        mpz_class temp = _y1 + _y2;
        mpz_mod(temp.get_mpz_t(), temp.get_mpz_t(), _p.get_mpz_t());

        if (temp == 0) {
            return point(0, 1, 0);
        }
        mpz_class first = 3 * _x1 * _x1 + _a;
        mpz_mod(first.get_mpz_t(), first.get_mpz_t(), _p.get_mpz_t()); // first = 3*x1^2 + a mod(_p)

        mpz_class second = 2 * _y1;
        mpz_invert(second.get_mpz_t(), second.get_mpz_t(), _p.get_mpz_t()); // second = (2y1)-^-1 mod(_p)

        m = first * second;
        mpz_mod(m.get_mpz_t(), m.get_mpz_t(), _p.get_mpz_t());

    }
    else {
        mpz_class first = _y2 - _y1;  // first = y2 - y1
        mpz_class second = _x2 - _x1; 
           
        mpz_invert(second.get_mpz_t(), second.get_mpz_t(), _p.get_mpz_t()); // second = (x2 - x1)^-1

        m = first * second;
        mpz_mod(m.get_mpz_t(), m.get_mpz_t(), _p.get_mpz_t());

    }
    mpz_class res_x = m * m - _x1 - _x2;
    mpz_mod(res_x.get_mpz_t(), res_x.get_mpz_t(), _p.get_mpz_t()); // x3 = m^2 - x1 -x2 mod(_p) 

    mpz_class res_y = m * (_x1 - res_x) - _y1; // y3 = m * (x1 - x3) - y1 mod(_p)
    mpz_mod(res_y.get_mpz_t(), res_y.get_mpz_t(), _p.get_mpz_t());

    return point(res_x, res_y, 1);;
}

std::ostream& operator<<(std::ostream& out, const elliptic_curve::_Point& P)  {
    if (P.get_z() == false) {
        out << "O\n";
    }
    else {
        out << "( " << P.get_x() << ", " << P.get_y() << " )";
    }
    return out;
}

elliptic_curve::_Point elliptic_curve::x2(const elliptic_curve::_Point& P) const {
    return sum(P, P);
}

elliptic_curve::_Point elliptic_curve::mul(const elliptic_curve::_Point& P, const mpz_class n) const {
    if (!P.get_z()) {
        return P;
    }
    mpz_class abs_n = abs(n);
    bool negative = false;
    if (n == 0) { // 0 * P = 0
        return new_point(0, 1, false);
    }
    else if (n < 0) { // -n*P = (|n|P)^-1        
        negative = true;
    }
    elliptic_curve::_Point result = new_point(0, 1,false);
    elliptic_curve::_Point temp = new_point(P.get_x(), P.get_y());
    auto bits = utils::binary(abs_n); // inverse bits of n
    for (auto& bit : bits) {
        if (bit == '1') {
            result = sum(result, temp); // result = result + temp
        }
        temp = x2(temp); // temp = temp + temp
    }
    if (negative) {
        temp = result;
        result = neg(temp); // (|n|P)^-1
    }

    mpz_mod(result.get_y().get_mpz_t(), result.get_y().get_mpz_t(), _p.get_mpz_t());
    return result;
}

elliptic_curve::_Point elliptic_curve::mul(const elliptic_curve::_Point& P, const long int n) const {
    return mul(P, mpz_class(n));
}

elliptic_curve::_Point elliptic_curve::neg(const elliptic_curve::_Point& P) const  {
    if (P.get_y() != 0) {
        elliptic_curve::_Point result = new_point(P.get_x(), -P.get_y());
        mpz_mod(result.get_y().get_mpz_t(), result.get_y().get_mpz_t(), _p.get_mpz_t());
        return result;
    }
    else {
        return P;
    }
}

elliptic_curve::_Point elliptic_curve::generate_point() {

    gmp_randstate_t rand_state;
    gmp_randinit_mt(rand_state);

    mpz_class x, t;
    gmp_randseed_ui(rand_state, utils::get_seed());
    
    while (true) {
        mpz_urandomm(x.get_mpz_t(), rand_state, _p.get_mpz_t());
        t = x * x * x + _a * x + _b;
        mpz_mod(t.get_mpz_t(), t.get_mpz_t(), _p.get_mpz_t()); //t = (x(x^2 + a) + b (mod p)

        if (mpz_legendre(t.get_mpz_t(), _p.get_mpz_t()) == -1) {
            continue;
        }
        else
            break;
    }
    t = utils::sqrtm(t, _p); // t = sqrt(t) in mod (_p)
    return point(x, t);
}

elliptic_curve::_Point elliptic_curve::generate_point(mpz_class x) {
    mpz_mod(x.get_mpz_t(), x.get_mpz_t(), _p.get_mpz_t());
    
    mpz_class t = x * x * x + _a * x + _b;
    mpz_mod(t.get_mpz_t(), t.get_mpz_t(), _p.get_mpz_t());


    if (mpz_legendre(t.get_mpz_t(), _p.get_mpz_t()) == -1) {
        return point (0, 1, false);
    }

    t = utils::sqrtm(t, _p); // t = sqrt(t) in mod (_p)
    return point(x, t);
}

std::vector<mpz_class> cross(std::vector<mpz_class> A, std::vector <mpz_class> B) {
    std::sort(A.begin(), A.end(), [](mpz_class first, mpz_class second) {
        if (mpz_cmp(first.get_mpz_t(), second.get_mpz_t()) < 0)
            return true;
        return false;
        });

    std::sort(B.begin(), B.end());


    unsigned long i = 0, j = 0;
    std::vector<mpz_class> S;

    while ((i < A.size()) && (j < B.size())) {
        if (A[i] <= B[j]) {
            if (A[i] == B[j]) {
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
    return S;

}

std::vector<mpz_class> shanks(const elliptic_curve& E,
    std::vector<mpz_class>& A, std::vector<mpz_class>& B, 
    const elliptic_curve::_Point P, const mpz_class W, const mpz_class x
) {
    mpz_class idx = 0, bound = W - 1;
    A.reserve(bound.get_ui());
    elliptic_curve::_Point point_temp = E.new_point(0, 1 , false);
   
    for (; idx <= bound; idx += 1) {
        point_temp = E.mul(P, E.get_p() + idx + 1);
        A.push_back(point_temp.get_x());
    }

    idx = 0;
    bound = W;
    A.reserve(bound.get_ui());
    
    for (; idx <= bound; idx += 1) {
        point_temp = E.mul(P, W * idx);
        B.push_back(point_temp.get_x());
    }
    return cross(A, B);
}

mpz_class ind(std::vector<mpz_class>& S, mpz_class s) {
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

mpz_class return_t(const elliptic_curve::_Point& P, elliptic_curve E, const mpz_class& beta, const mpz_class& gamma, const mpz_class& W, const mpz_class p) {
    mpz_class temp = beta + gamma * W;
    mpz_class test = p + 1 + temp;
    auto pnt = E.mul(P, test);
    if (pnt.get_z() == false)
        return temp;
    else {
        return beta - gamma * W;
    }
}

mpz_class order(const elliptic_curve& curve) {
    if (curve.get_p() <= 229) {

    }


    elliptic_curve E(curve);
    mpz_class g = utils::rand_not_sqr_res(curve.get_p());

    mpz_class W;
    mpf_class temp = curve.get_p();
    temp = sqrt(curve.get_p());
    temp = sqrt(temp);
    temp *= sqrt(2);
    mpf_ceil(temp.get_mpf_t(), temp.get_mpf_t());
    W = temp;


    mpz_class c = g * g * curve.get_a();
    mpz_class d = g * g * g * curve.get_b();

    mpz_class sigma, y_2, a_x;
    gmp_randstate_t rand_state;
    gmp_randinit_mt(rand_state);

    while (true) {

        mpz_class x;

        gmp_randseed_ui(rand_state, utils::get_seed());
        mpz_urandomm(x.get_mpz_t(), rand_state, curve.get_p().get_mpz_t());

        mpz_class temp = x * x * x + curve.get_a() * x + curve.get_b();
        mpz_mod(temp.get_mpz_t(), temp.get_mpz_t(), curve.get_p().get_mpz_t());
        short int sigma = mpz_legendre(temp.get_mpz_t(), curve.get_p().get_mpz_t());

        if (sigma == 0) {
            continue;
        }
        if (sigma == 1) {
            ;
        }
        else {

            E.set_a(c);
            E.set_b(d);
            x *= g;

        }
        auto P = E.generate_point(x);
        if (P.get_z() == false) {
            E.set_a(curve.get_a());
            E.set_b(curve.get_b());
            continue;
        }

        std::vector<mpz_class> A, B, S;

        S = shanks(E, A, B, P, W, x);

        if (S.size() != 1) {
            E.set_a(curve.get_a());
            E.set_b(curve.get_b());
            continue;
        }

        gmp_randclear(rand_state);
        auto s = S[0];
        auto beta = ind(A, s);
        auto gamma = ind(B, s);

        auto t = return_t(P, E, beta, gamma, mpz_class(W), mpz_class(E.get_p()));
        return mpz_class(E.get_p() + 1 + sigma * t);
    }

}

int main() {

    mpz_class a = 2, b = 3, p("421");//p("599234844171323798014378139219684288087362388635279936453684278910360928027555754598727507571581246001777277670122448689713029876360082601445563780429103017907332764522609770882588110158952861810496271751306659999739053379195854324818229379783745932907778328817332573106903755814283051471753305325707");
    mpz_class c = 120, d = 1, e = 327, f = 12421;


    elliptic_curve curve(a, b, p);
   /* auto pnt1 = curve.generate_point();
    auto pnt2 = curve.generate_point();
    for (int i = 0; i < 100000; ++i) {

        //auto pnt3 = curve.sum(pnt1, pnt2);
        auto pnt4 = curve.mul(pnt1, 531432);
        //pnt3 = curve.neg(pnt3);
    }*/
    std::cout << order(curve);
/*
    point pnt1 = curve.new_point(3, 6);
    point pnt2 = curve.new_point(80, -10);
*/
    //point pnt = curve.generate_point();
    //std::cout << pnt << std::endl;


}



//mpz_set_str(p, "599234844171323798014378139219684288087362388635279936453684278910360928027555754598727507571581246001777277670122448689713029876360082601445563780429103017907332764522609770882588110158952861810496271751306659999739053379195854324818229379783745932907778328817332573106903755814283051471753305325707", 10);

