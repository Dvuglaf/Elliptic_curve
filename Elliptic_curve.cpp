#include <iostream>
#include "Elliptic_curve.h"
#include "utils.h"
#include <algorithm>
#include <gmpxx.h>

mpz_class elliptic_curve::_Point::get_x() const {
    return _x;
}

mpz_class elliptic_curve::_Point::_Point::get_y() const {
    return _y;
}

bool elliptic_curve::_Point::_Point::get_z() const {
    return _z;
}

elliptic_curve::_Point::_Point(const mpz_class& x, const mpz_class& y, const bool z = true) : _x(x), _y(y), _z(z) {}

elliptic_curve::elliptic_curve(const mpz_class& a, const mpz_class& b, const mpz_class& p) : _a(a), _b(b), _p(p) {
    if (-16 * (4 * a * a * a + 27 * _b * _b) != 0) {
        ;
    }
    else {
        std::cout << "Discriminant equals zero. Setting a = 0, b = 1, p = 5.\n";
        _a = 0;
        _b = 1;
        _p = 5;
    }

    if (_p < 5) {
        std::cout << "Low module. Setting p = 5.\n";
    }
}

void elliptic_curve::set_a(const mpz_class& a) {
    if (-16 * (4 * a * a * a + 27 * _b * _b) != 0) {
        _a = a;
    }
    else {
        std::cout << "Discriminant equals zero. No changes.\n";
    }
}
void elliptic_curve::set_b(const mpz_class& b) {
    if (-16 * (4 * _a * _a * _a + 27 * b * b) != 0) {
        _b = b;
    }
    else {
        std::cout << "Discriminant equals zero\n";
    }
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

bool elliptic_curve::exist_point(const mpz_class& x, const mpz_class& y, const bool z = true) const {
    if (!z)
        return true;

    mpz_class left = y * y;
    mpz_mod(left.get_mpz_t(), left.get_mpz_t(), _p.get_mpz_t());
    mpz_class right = x * x * x + x * _a + _b;
    mpz_mod(right.get_mpz_t(), right.get_mpz_t(), _p.get_mpz_t());

    if (left != right) {
        std::cout << "Point = (" << x << ", " << y << ") are not on the curve.\n";
        return false;
    }
    return true;
}

elliptic_curve::_Point elliptic_curve::new_point(const mpz_class& x, const mpz_class& y, const bool z = true) const {
    if (!z) {
        return point(0, 1, false);
    }
    if (elliptic_curve::exist_point(x, y, z)) {
        auto copy_x = x;
        auto copy_y = y;
        mpz_mod(copy_x.get_mpz_t(), copy_x.get_mpz_t(), _p.get_mpz_t());
        mpz_mod(copy_y.get_mpz_t(), copy_y.get_mpz_t(), _p.get_mpz_t());

        return point(copy_x, copy_y, z);
    }
    return point(0, 1, false);
}

elliptic_curve::_Point elliptic_curve::new_point(const long int x, const long int y, const bool z = true) const {
    return new_point(mpz_class(x), mpz_class(y), z);
}

elliptic_curve::_Point elliptic_curve::sum(const elliptic_curve::_Point& P, const elliptic_curve::_Point& Q) const {
    const auto _x1 = P.get_x();
    const auto _y1 = P.get_y();
    const auto _x2 = Q.get_x();
    const auto _y2 = Q.get_y();

    if (P.get_z() == 0) return Q; // P = O
    if (Q.get_z() == 0) return P; // Q = O

    mpz_class m;

    if (_x1 == _x2) {
        mpz_class y_sum = _y1 + _y2;
        mpz_mod(y_sum.get_mpz_t(), y_sum.get_mpz_t(), _p.get_mpz_t());

        if (y_sum == 0) {
            return point(0, 1, false);
        }
        mpz_class first = 3 * _x1 * _x1 + _a;
        mpz_mod(first.get_mpz_t(), first.get_mpz_t(), _p.get_mpz_t()); // first = 3*x1^2 + a mod(_p)

        mpz_class second = 2 * _y1;
        mpz_invert(second.get_mpz_t(), second.get_mpz_t(), _p.get_mpz_t()); // second = (2y1)-^-1 mod(_p)

        m = (first * second) % _p;
    }
    else {
        const mpz_class first = _y2 - _y1;  // first = y2 - y1
        mpz_class second = _x2 - _x1; 
           
        mpz_invert(second.get_mpz_t(), second.get_mpz_t(), _p.get_mpz_t()); // second = (x2 - x1)^-1 mod(_p)

        m = (first * second) % _p;

    }
    mpz_class res_x = m * m - _x1 - _x2;
    mpz_mod(res_x.get_mpz_t(), res_x.get_mpz_t(), _p.get_mpz_t()); // x3 = m^2 - x1 -x2 mod(_p) 

    mpz_class res_y = m * (_x1 - res_x) - _y1; 
    mpz_mod(res_y.get_mpz_t(), res_y.get_mpz_t(), _p.get_mpz_t()); // y3 = m * (x1 - x3) - y1 mod(_p)

    return point(res_x, res_y, 1);
}

elliptic_curve::_Point elliptic_curve::sub(const elliptic_curve::_Point& P, const elliptic_curve::_Point& Q) const {
    return sum(P, neg(Q));
}

std::ostream& operator<<(std::ostream& out, const elliptic_curve::_Point& P)  {
    if (P.get_z() == false) {
        out << "O";
    }
    else {
        out << "( " << P.get_x() << ", " << P.get_y() << " )";
    }
    return out;
}

elliptic_curve::_Point elliptic_curve::x2(const elliptic_curve::_Point& P) const {
    return sum(P, P);
}

elliptic_curve::_Point elliptic_curve::mul1(const elliptic_curve::_Point& P, const mpz_class& n) const {
    if (!P.get_z()) {
        return P;
    }
    const mpz_class abs_n = abs(n);
    bool negative = false;
    if (n == 0) { // 0 * P = 0
        return new_point(0, 1, false);
    }
    else if (n < 0) { // -n*P = (|n|P)^-1        
        negative = true;
    }
    auto result = new_point(0, 1, false);
    auto temp = new_point(P.get_x(), P.get_y()); // temp = P

    const auto bits = utils::binary(abs_n); // inverse bits of n
    
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

    mpz_mod(result._x.get_mpz_t(), result._x.get_mpz_t(), _p.get_mpz_t());
    mpz_mod(result._y.get_mpz_t(), result._y.get_mpz_t(), _p.get_mpz_t());

    return result;
}

elliptic_curve::_Point elliptic_curve::mul2(const elliptic_curve::_Point& P, const mpz_class& n) const {
    if (n == 0) {
        return _Point(0, 1, false);
    }
    auto result = P;
    const auto bits_m = utils::binary(3 * n);
    const auto bits_n = utils::binary_and_n_zeros_add(n, bits_m.size());

    for (unsigned long j = bits_m.size() - 2; j >= 1; --j) {
        result = x2(result);
        if ((bits_m[j] == '1') && (bits_n[j] == '0')) {
            result = sum(result, P);
        }
        if ((bits_m[j] == '0') && (bits_n[j] == '1')) {
            result = sub(result, P);
        }
    }

    mpz_mod(result._x.get_mpz_t(), result._x.get_mpz_t(), _p.get_mpz_t());
    mpz_mod(result._y.get_mpz_t(), result._y.get_mpz_t(), _p.get_mpz_t());

    return result;
}

elliptic_curve::_Point elliptic_curve::neg(const elliptic_curve::_Point& P) const  {
    if (P.get_y() != 0) {
        auto result = new_point(P.get_x(), -P.get_y());
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
    
    mpz_class t = (x * x * x + _a * x + _b) % _p;

    if (mpz_legendre(t.get_mpz_t(), _p.get_mpz_t()) == -1) {
        return point (0, 1, false);
    }

    t = utils::sqrtm(t, _p); // t = sqrt(t) in mod (_p)
    return point(x, t);
}

std::vector<mpz_class> cross(std::vector<mpz_class> A, std::vector <mpz_class> B) {
    std::sort(A.begin(), A.end());

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
                S.resize(0);
                return S;
            }
        }
        else {
            ++j;
            while ((j < B.size() - 1) && (B[j] == B[j - 1])) {
                ++j;
                S.resize(0);
                return S;
            }
        }
    }
    return S;
}

std::vector<mpz_class> shanks(
    const elliptic_curve& E, std::vector<mpz_class>& A, std::vector<mpz_class>& B, 
    const elliptic_curve::_Point& P, const mpz_class& W) {

    mpz_class idx = 0, bound = W - 1;
    auto point_temp = E.mul1(P, E.get_p());

    for (; idx <= bound; idx += 1) {
        point_temp = E.sum(point_temp, P);
        A.push_back(point_temp.get_x());
    }

    idx = 0;
    bound = W;
    B.push_back(0);
    const auto pnt = E.mul1(P, W);
    B.push_back(pnt.get_x());
    for (idx = 2; idx <= bound; idx += 1) {
        point_temp = E.sum(point_temp, pnt);
        B.push_back(point_temp.get_x());
    }
    return cross(A, B);
}

mpz_class ind(const std::vector<mpz_class>& S, const mpz_class& s) {
    mpz_class i = 0;
    for (auto& it : S) {
        if (it == s) {
            return i;
        }
        i += 1;
    }
    return mpz_class(-1);

}

mpz_class return_t(const elliptic_curve::_Point& P, const elliptic_curve& E, const mpz_class& beta, const mpz_class& gamma, const mpz_class& W, const mpz_class& p) {
    const mpz_class temp = beta + gamma * W;
    const mpz_class test = p + 1 + temp;
    auto pnt = E.mul1(P, test);
    if (pnt.get_z() == false)
        return temp;
    else {
        return beta - gamma * W;
    }
}

mpz_class order(const elliptic_curve& curve) {
    if (curve.get_p() <= 229) {
        mpz_class nominator, result = curve.get_p() + 1;
        for (int x = 0; x < curve.get_p(); ++x) {
            nominator = x * x * x + curve.get_a() * x + curve.get_b();
            mpz_mod(nominator.get_mpz_t(), nominator.get_mpz_t(), curve.get_p().get_mpz_t());
            result += mpz_legendre(nominator.get_mpz_t(), curve.get_p().get_mpz_t());
        }
        return result;
    }

    elliptic_curve E(curve);
    const mpz_class g = utils::rand_not_sqr_res(curve.get_p());

    mpz_class W;
    mpf_class temp = curve.get_p();
    temp = sqrt(curve.get_p());
    temp = sqrt(temp);
    temp *= sqrt(2);
    mpf_ceil(temp.get_mpf_t(), temp.get_mpf_t());
    W = temp;


    const mpz_class c = g * g * curve.get_a();
    const mpz_class d = g * g * g * curve.get_b();

    gmp_randstate_t rand_state;
    gmp_randinit_mt(rand_state);

    while (true) {

        mpz_class x;

        gmp_randseed_ui(rand_state, utils::get_seed());
        mpz_urandomm(x.get_mpz_t(), rand_state, curve.get_p().get_mpz_t());

        mpz_class nominator = x * x * x + curve.get_a() * x + curve.get_b();
        mpz_mod(nominator.get_mpz_t(), nominator.get_mpz_t(), curve.get_p().get_mpz_t());
        
        const int sigma = mpz_legendre(nominator.get_mpz_t(), curve.get_p().get_mpz_t());

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
        
        const auto P = E.generate_point(x);
        if (P.get_z() == false) {
            E.set_a(curve.get_a());
            E.set_b(curve.get_b());
            continue;
        }

        std::vector<mpz_class> A, B;

        const auto S = shanks(E, A, B, P, W);

        if (S.size() != 1) {
            E.set_a(curve.get_a());
            E.set_b(curve.get_b());
            continue;
        }

        gmp_randclear(rand_state);
        const auto s = S[0];
        const auto beta = ind(A, s);
        const auto gamma = ind(B, s);

        const auto t = return_t(P, E, beta, gamma, W, E.get_p());
        return mpz_class(E.get_p() + 1 + sigma * t);
    }
}

int main() {

    mpz_class a = 2, b = 3, p("599234844171323798014378139219684288087362388635279936453684278910360928027555754598727507571581246001777277670122448689713029876360082601445563780429103017907332764522609770882588110158952861810496271751306659999739053379195854324818229379783745932907778328817332573106903755814283051471753305325707");//p("599234844171323798014378139219684288087362388635279936453684278910360928027555754598727507571581246001777277670122448689713029876360082601445563780429103017907332764522609770882588110158952861810496271751306659999739053379195854324818229379783745932907778328817332573106903755814283051471753305325707");

    elliptic_curve curve(a, b, p);
    //std::cout << curve.generate_point();

    auto pnt1 = curve.new_point(mpz_class("289220140802772654815015008848224377106589756640499318037321215985542245573510325194905346787511430150282038865361627058777058918115664720836662973489139462074770880145107943861816555494608423600911860523394797999291461136771320750592578542164638864581223568762917722592396298086270572433888180389227"),
                                mpz_class("283139674163665446610040439261624655768751271240266493137719217364391736830777551106167054698249968822785660388608546719559141814330389549252902429782513111562645687295347496315545318709925367966950167844732317643888279130540366236692918982910769248252048230954794430247412224105614992699080004620485"));

    auto pnt2 = curve.new_point(mpz_class("278681615526238039257096200453748761033454289735737304845430341131529189020613844951822326569210198117364983609734402420521182281026677393074644225105217587258056132553670682034462158607098943830531565904798784646782152450387011329669045620962567222200097053774069826512866415100137946624421005304236"),
                                mpz_class("70382567546435466398715735681811268162879856507961739298500901726423490122795241253090483191058744185447902911288049620355364708699715879988728683722865480567759635975627343148410254584765121659806097096082477610296996473605244245778834021419434521481982215713580628333489076004466991353163461651730"));
    
    std::cout << curve.sum(pnt1, pnt2) << std::endl;
    std::cout << curve.sub(pnt1, pnt2) << std::endl;
    std::cout << curve.mul1(pnt1, 412) << std::endl;
    std::cout << curve.mul2(pnt2, mpz_class("194120")) << std::endl << std::endl;
    
    elliptic_curve curve_test_low(a, b, 421);
    pnt1 = curve_test_low.new_point(308, 85);
    pnt2 = curve_test_low.new_point(380, 329);
    std::cout << curve_test_low.sum(pnt1, pnt2) << std::endl;
    std::cout << curve_test_low.sub(pnt1, pnt2) << std::endl;
    std::cout << curve_test_low.mul1(pnt1, 41) << "\t\t";
    std::cout << curve_test_low.mul2(pnt1, 41) << std::endl;
    std::cout << curve_test_low.mul1(pnt2, 41) << "\t\t";
    std::cout << curve_test_low.mul2(pnt2, 41) << std::endl;

    std::cout << order(curve_test_low);


    /*
E = EllipticCurve(
        GF(599234844171323798014378139219684288087362388635279936453684278910360928027555754598727507571581246001777277670122448689713029876360082601445563780429103017907332764522609770882588110158952861810496271751306659999739053379195854324818229379783745932907778328817332573106903755814283051471753305325707),
        [2,3]
    );
P = E(289220140802772654815015008848224377106589756640499318037321215985542245573510325194905346787511430150282038865361627058777058918115664720836662973489139462074770880145107943861816555494608423600911860523394797999291461136771320750592578542164638864581223568762917722592396298086270572433888180389227, 
      283139674163665446610040439261624655768751271240266493137719217364391736830777551106167054698249968822785660388608546719559141814330389549252902429782513111562645687295347496315545318709925367966950167844732317643888279130540366236692918982910769248252048230954794430247412224105614992699080004620485
     );

Q = E(278681615526238039257096200453748761033454289735737304845430341131529189020613844951822326569210198117364983609734402420521182281026677393074644225105217587258056132553670682034462158607098943830531565904798784646782152450387011329669045620962567222200097053774069826512866415100137946624421005304236, 
      70382567546435466398715735681811268162879856507961739298500901726423490122795241253090483191058744185447902911288049620355364708699715879988728683722865480567759635975627343148410254584765121659806097096082477610296996473605244245778834021419434521481982215713580628333489076004466991353163461651730
     );

print(P + Q);
print(P - Q);
print(P * 412);
print(Q * 194120);
    
    */
    /*
E = EllipticCurve(GF(421), [2, 3]);
P = E(308, 85);
Q = E(380, 329);

print(P+Q);
print(P-Q);
print(41*P);
print(41*Q);
print(order(E));
    */
}



//mpz_set_str(p, "599234844171323798014378139219684288087362388635279936453684278910360928027555754598727507571581246001777277670122448689713029876360082601445563780429103017907332764522609770882588110158952861810496271751306659999739053379195854324818229379783745932907778328817332573106903755814283051471753305325707", 10);

