#pragma once
#include <gmp.h>

namespace utils {
	void sqrtm(mpz_t, const mpz_t, const mpz_t);

	std::string binary(const mpz_t);

	void rand_not_sqr_res(mpz_t, const mpz_t);

	unsigned long get_seed();
}