#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <fstream>

using namespace std;

int main()
{
	char inputStr[1024];

	// input a
	mpz_t a;
	cout << "Enter a: ";
	cin >> inputStr;
	mpz_init_set_str(a, inputStr, 10);

	// input c
	unsigned long c;
	cout << "Enter c: ";
	cin >> c;

	// input m2exp
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

	// open file
	FILE *fileOut;
	char fileNameOut[] = "data.txt";
	if ((fileOut = fopen(fileNameOut, "w")) == NULL)
	{
		cout << "error opening fileOut " << fileNameOut << endl;
		return 1;
	}

	// input ArraySize
	unsigned long int ArraySize;
	cout << "Enter ArraySize: ";
	cin >> ArraySize;

	// Random Number Generation [see 5.13]
	mpz_t x; // stores the generated numbers
	mpz_init(x);
	for (unsigned long int i = 0; i < ArraySize; i++)
	{
		mpz_urandomb(x, state, m2exp); // generates next random number
		mpz_out_str(fileOut, 10, x);   // write to file
		fprintf(fileOut, "\n");		   // separated by '\n'
	}
	fclose(fileOut);
	cout << "---the end---" << endl;

	return 0;
}
