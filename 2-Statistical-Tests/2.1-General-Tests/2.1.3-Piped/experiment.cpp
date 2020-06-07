#include <iostream>
#include <cstdio>
#include <gmp.h>
#include <unistd.h>

using namespace std;

int main()
{
    char inputStr[1024];

	// input a: multiplier
	mpz_t a;
	cout << "Enter a: ";
	cin >> inputStr;
	mpz_init_set_str(a, inputStr, 10);

	// input c: increment
	unsigned long c;
	cout << "Enter c: ";
	cin >> c;

	// input m2exp: modulus
	mp_bitcnt_t m2exp;
	cout << "Enter m2exp: ";
	cin >> m2exp;

	// Random State Initialization [see 9.1]
	gmp_randstate_t state;
	gmp_randinit_lc_2exp(state, a, c, m2exp);

	// Random State Seeding [see 9.2]
	mpz_t seed;
	cout << "Enter seed: ";
	cin >> inputStr;
	mpz_init_set_str(seed, inputStr, 10);
	gmp_randseed(state, seed);

	mpz_t xZ; // holds the generated RNs
	mpz_init(xZ);

	if (pipe(p) < 0)
	{
		cout << "Error in piping!" << endl;
		exit(1);
	}
	
	for (unsigned long int i = 0; i < ArraySize; i++)
	{
		mpz_urandomb(xZ, state, m2exp);
		mpz_get_str(inputStr, 10, xZ);
		write(p[1], inputStr, 1024);
		read(p[0], inputStr, 1024);
		cout<<inputStr<<endl;
	}
    
    return 0;
}