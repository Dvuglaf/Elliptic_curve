#include <random>
#include "utils.h"

void utils::rand_not_sqr_res(mpz_t result, const mpz_t p) {
    mpz_t temp;
    mpz_init_set(temp, p);
    mpz_cdiv_q_ui(temp, temp, 2); // temp = p/2
    while (true) {
        if (mpz_legendre(temp, p) == -1) {
            mpz_set(result, temp);
            mpz_clear(temp);
            return;
        }
        mpz_add_ui(temp, temp, 1);
    }

}

unsigned long utils::get_seed() {
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_int_distribution<unsigned long> uniform_dist(1, std::numeric_limits<unsigned long>::max());
    unsigned long value = uniform_dist(e1);
    return value;

}

std::string utils::binary(const mpz_t x) {
    auto len = mpz_sizeinbase(x, 2) + 2;
    char* c_str = new char[len];
    c_str = mpz_get_str(c_str, 2, x);
    std::string str(c_str);
    delete[] c_str;
    std::reverse(str.begin(), str.end());
    return str;
}

void utils::sqrtm(mpz_t result, const mpz_t x, const mpz_t mod) {
    mpz_t a, p;
    mpz_init_set(a, x);
    mpz_init_set_ui(p, 8);
    mpz_mod(a, a, mod);

    mpz_mod_ui(p, mod, 8); // p = _p mod 8

    // 1. [Простейшие случаи]
    if (mpz_cmp_ui(p, 3) == 0 || mpz_cmp_ui(p, 7) == 0) {
        mpz_t pow;
        mpz_init_set(pow, mod);
        mpz_add_ui(pow, pow, 1); // pow = _p + 1
        mpz_cdiv_q_ui(pow, pow, 4); // pow = (_p + 1)/4
        mpz_powm(result, a, pow, mod); // result = a^pow mod(_p)
        mpz_clears(a, p, pow, NULL);
        return;
    }
    if (mpz_cmp_ui(p, 5) == 0) {
        mpz_t pow, c;
        mpz_init_set(pow, mod);
        mpz_add_ui(pow, pow, 3);
        mpz_cdiv_q_ui(pow, pow, 8); // pow = (_p + 3)/8
        mpz_powm(result, a, pow, mod); // result = a^pow mod(_p)
        mpz_init(c);
        mpz_powm_ui(c, result, 2, mod); // c = x^2 mod(_p)
        if (mpz_cmp(c, a) != 0) {
            mpz_set(pow, mod);
            mpz_sub_ui(pow, pow, 1);
            mpz_cdiv_q_ui(pow, pow, 4); // pow = (_p - 1)/4
            mpz_t second;
            mpz_init_set_ui(second, 2);
            mpz_powm(second, second, pow, mod); // second = 2^pow mod(_p)
            mpz_mul(result, result, second);
            mpz_mod(result, result, mod); // result = result * second mod(_p)
            mpz_clears(second, c, a, p, pow, NULL);
            return;
        }
        mpz_clears(c, a, p, pow, NULL);
        return;
    }

    // 2. [Случай _p = 1 (mod 8)]
    //находим случ. квадратичный невычет (legandre = - 1)

    mpz_t d, t, q;
    mpz_inits(d, t, q, NULL);
    rand_not_sqr_res(d, mod); // d = random not square residue mod(_p)

    //present q = p-1 = t * 2^s
    mpz_set(q, mod);
    mpz_sub_ui(q, q, 1);
    unsigned int s = 0;
    while (mpz_tstbit(q, s) == 0) s++;
    mpz_t second;
    mpz_init_set_ui(second, 2);
    mpz_pow_ui(second, second, s); // second = 2^s
    mpz_div(t, q, second); // t = (p-1)/2^s
    mpz_clears(q, p, NULL);


    mpz_t A, D, m;
    mpz_inits(A, D, m, NULL);

    mpz_powm(A, a, t, mod); // A = a^t mod(_p)
    mpz_powm(D, d, t, mod); // D = d^t mod(_p)

    mpz_set_ui(m, 0);
    mpz_t base, pow1, pow2;
    mpz_inits(base, pow1, pow2, NULL);
    for (unsigned int i = 0; i < s; ++i) {
        mpz_set(base, D);
        mpz_powm(base, base, m, mod);
        mpz_mul(base, A, base);
        mpz_mod(base, base, mod); // base = A * D^m mod(_p)
        mpz_set_ui(pow1, 2);
        mpz_set_ui(pow2, s);
        mpz_sub_ui(pow2, pow2, 1);
        mpz_sub_ui(pow2, pow2, i); // pow2 = s - 1 - i

        mpz_powm(pow1, pow1, pow2, mod); // pow1 = 2^pow2 mod(_p)
        mpz_powm(base, base, pow1, mod); // base = base^pow1 mod(_p)

        mpz_set(second, mod);
        mpz_sub_ui(second, second, 1); // second = -1 mod(_p)


        if (mpz_cmp(base, second) == 0) {
            mpz_set_ui(base, 2);
            mpz_set_ui(pow1, i);
            mpz_powm(base, base, pow1, mod);
            mpz_add(m, m, base); // m = m + 2^i
        }
    }
    mpz_set(pow1, t);
    mpz_add_ui(pow1, pow1, 1);
    mpz_div_ui(pow1, pow1, 2); // pow1 = (t+1)/2

    mpz_set(pow2, m);
    mpz_div_ui(pow2, pow2, 2); // pow2 = m/2 

    mpz_powm(a, a, pow1, mod); // a = a^pow1 mod(_p)
    mpz_powm(D, D, pow2, mod); // D = D^pow2 mod(_p)

    mpz_mul(result, a, D);
    mpz_mod(result, result, mod); //result = a * D mod(_p)
    mpz_clears(a, d, t, second, A, D, m, base, pow1, pow2, NULL);

}
