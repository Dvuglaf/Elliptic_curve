#include <random>
#include "utils.h"
#include <iostream>

mpz_class utils::rand_not_sqr_res(const mpz_class& p) {
    mpz_class result = p;
    result /= 2; // result = p/2
    while (true) {
        if (mpz_legendre(result.get_mpz_t(), p.get_mpz_t()) == -1) {
            return result;
        }
        result += 1;
    }
}

unsigned long utils::get_seed() {
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_int_distribution<unsigned long> uniform_dist(1, std::numeric_limits<unsigned long>::max());
    unsigned long value = uniform_dist(e1);
    return value;

}

std::string utils::binary(const mpz_class& x) {
    std::string str = x.get_str(2);
    std::reverse(str.begin(), str.end());
    return str;
}

std::string utils::ternary(const mpz_class& n) {
    if (n == 1)
        return std::string("1");
    std::string bits = utils::binary(n);
    for (int i = 0; i < bits.size() - 1; ++i) {
        if (bits[i] == '1' && bits[i + 1] == '1') {
            int t = 0;
            int s = i;
            while (bits[i] == '1' && i < bits.size()) {
                ++i;
                ++t;
            }
            bits.replace(s, 1, "A");
            bits.replace(s + 1, t - 1, std::string(t - 1, '0'));
            if (s + t < bits.size()) {
                bits.replace(s + t, 1, "1");
            }
            else {
                bits.insert(s + t, 1, '1');
            }
            --i;
        }

    }
    return bits;

}

std::string utils::binary_and_n_zeros_add(const mpz_class& x, const unsigned long n) {
    std::string str = x.get_str(2);
    std::reverse(str.begin(), str.end());
    while (str.size() < n) {
        str.push_back('0');
    }
    return str;
}

mpz_class utils::sqrtm(const mpz_class& x, const mpz_class& mod) {
    // Alg 2.3.8 

    mpz_class a = x, p = 8, result;

    mpz_mod(a.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());
    mpz_mod_ui(p.get_mpz_t(), mod.get_mpz_t(), 8); // p = _p mod 8

    // 1. [Simplest cases]
    if (p == 3 || p == 7) {
        mpz_class pow = (mod + 1) / 4;
        mpz_powm(result.get_mpz_t(), a.get_mpz_t(), pow.get_mpz_t(), mod.get_mpz_t()); // result = a^pow mod(_p)
        return result;
    }
    if (p == 5) {
        mpz_class pow = (mod + 3) / 8, c;

        mpz_powm(result.get_mpz_t(), a.get_mpz_t(), pow.get_mpz_t(), mod.get_mpz_t()); // result = a^pow mod(_p)

        mpz_powm_ui(c.get_mpz_t(), result.get_mpz_t(), 2, mod.get_mpz_t()); // c = x^2 mod(_p)
        if (c != a) {
            pow = (mod - 1) / 4;
            mpz_class base = 2;
            mpz_powm(base.get_mpz_t(), base.get_mpz_t(), pow.get_mpz_t(), mod.get_mpz_t()); // base = 2^pow mod(_p)
            result *= base;
            mpz_mod(result.get_mpz_t(), result.get_mpz_t(), mod.get_mpz_t()); // result = result * base mod(_p)
        }
        return result;
    }

    // 2. [Case _p = 1 (mod 8)]
    //find random not square residue (legandre = - 1)

    mpz_class d;
    d = rand_not_sqr_res(mod); // d = random not square residue mod(_p)

    //view q = p-1 = t * 2^s
    mpz_class q = mod - 1;

    unsigned int s = 0;
    while (mpz_tstbit(q.get_mpz_t(), s) == 0) ++s;

    mpz_class t = q, denominator = 2;

    mpz_pow_ui(denominator.get_mpz_t(), denominator.get_mpz_t(), s); // second = 2^s
    t /= denominator;

    mpz_class A, D, m = 0;

    mpz_powm(A.get_mpz_t(), a.get_mpz_t(), t.get_mpz_t(), mod.get_mpz_t()); // A = a^t mod(_p)
    mpz_powm(D.get_mpz_t(), d.get_mpz_t(), t.get_mpz_t(), mod.get_mpz_t()); // D = d^t mod(_p)

    mpz_class left, pow1, pow2;
    const mpz_class right = mod - 1;
    for (mpz_class i = 0; i < s; i += 1) {
        left = D;
        mpz_powm(left.get_mpz_t(), left.get_mpz_t(), m.get_mpz_t(), mod.get_mpz_t());
        left *= A;
        mpz_mod(left.get_mpz_t(), left.get_mpz_t(), mod.get_mpz_t()); // left = A * D^m mod(_p)
        pow1 = 2;
        pow2 = s - 1 - i; // pow2 = s - 1 - i

        mpz_powm(pow1.get_mpz_t(), pow1.get_mpz_t(), pow2.get_mpz_t(), mod.get_mpz_t()); // pow1 = 2^pow2 mod(_p)
        mpz_powm(left.get_mpz_t(), left.get_mpz_t(), pow1.get_mpz_t(), mod.get_mpz_t()); // left = left^pow1 mod(_p)
       
        if (left == right) {
            mpz_class base = 2;
            mpz_powm(base.get_mpz_t(), base.get_mpz_t(), i.get_mpz_t(), mod.get_mpz_t());
            m += base;
        }
    }
    pow1 = (t + 1) / 2; // pow1 = (t+1)/2

    pow2 = m / 2; // pow2 = m/2 

    mpz_powm(a.get_mpz_t(), a. get_mpz_t(), pow1.get_mpz_t(), mod.get_mpz_t()); // a = a^pow1 mod(_p)
    mpz_powm(D.get_mpz_t(), D.get_mpz_t(), pow2.get_mpz_t(), mod.get_mpz_t()); // D = D^pow2 mod(_p)

    result = a * D;
    mpz_mod(result.get_mpz_t(), result.get_mpz_t(), mod.get_mpz_t()); //result = a * D mod(_p)
    return result;

}