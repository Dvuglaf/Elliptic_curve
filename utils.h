#pragma once
#include <gmpxx.h>

namespace utils {
	mpz_class sqrtm(const mpz_class&, const mpz_class&);

	std::string binary(const mpz_class&);
	
	std::string binary_and_n_zeros_add(const mpz_class&, const unsigned long);

	mpz_class rand_not_sqr_res(const mpz_class&);

	unsigned long get_seed();
}