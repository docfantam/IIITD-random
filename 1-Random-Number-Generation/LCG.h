#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

char* __rng_name_lcg()
{
	return "LCG";
}

void __rng_run_lcg(char* pname, char* fname, char* nname)
{
	if(pname)
	{
		mpz_t a;	// a
		mpz_init_set_str(a, strtok(pname, ":"), 10);

		unsigned long c;	// c
		c = atol(strtok(NULL, ":"));

		mp_bitcnt_t m2exp;	//m
		m2exp = atoi(strtok(NULL, ":"));

		// Random State Initialization [see 9.1]
		gmp_randstate_t state;
		gmp_randinit_lc_2exp(state, a, c, m2exp);

		// Random State Seeding [see 9.2]
		mpz_t seed;
		mpz_init_set_str(seed, strtok(NULL, ":"), 10);
		gmp_randseed(state, seed);

		// Random Number Generation [see 5.13]
		mpz_t x; // stores the generated numbers
		mpz_init(x);
		unsigned long ArraySize = 1000;
		if(nname)
			ArraySize = atol(nname);

		if(fname)	// FILE
		{
			// open file
			FILE *fileOut;
			char fileNameOut[100];
			strcpy(fileNameOut, fname);
			if ((fileOut = fopen(fileNameOut, "w")) == NULL)
			{
				printf("error opening fileOut %s\n", fileNameOut);
				return;
			}
			
			for (unsigned long int i = 0; i < ArraySize; i++)
			{
				mpz_urandomb(x, state, m2exp); // generates next random number
				mpz_out_str(fileOut, 10, x);   // write to file
				fprintf(fileOut, "\n");		   // separated by '\n'
			}
			fclose(fileOut);
		}
		else	// STDOUT
		{
			for (unsigned long int i = 0; i < ArraySize; i++)
			{
				mpz_urandomb(x, state, m2exp); // generates next random number
				gmp_printf("%Zd\n", x);   // write to file
			}
		}
	}
	else 
	{
		printf("parameter is a required field!\n");
	}
}

char* __rng_desc_lcg()
{
	return "Linear Congruential Generator\n\
	The random numbers are generated as X = (a*X + c) mod m where \n\
	a: multiplier\n\
	c: increment\n\
	m: modulus\n\n\
	Usage: ./a.out -g LCG -p <a>:<c>:<m>:<seed> [-n size] [-f filename]\n";
}