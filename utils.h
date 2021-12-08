#pragma once
#include <gmpxx.h>

namespace utils {
	mpz_class sqrtm(const mpz_class&, const mpz_class&);

	std::string binary(const mpz_class&);

	mpz_class rand_not_sqr_res(const mpz_class&);

	unsigned long get_seed();
}