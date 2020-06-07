#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <fstream>
#include <boost/math/special_functions/gamma.hpp> // for normalized incomplete gamma function

using namespace std;

int main()
{
	char inputStr[1024];
	unsigned long int i; // loop iterator

	// input m2exp [Note: the range of random numbers generated is 0 to 2^m2exp - 1]
	mp_bitcnt_t m2exp;
	cout << "Enter m2exp: ";
	cin >> m2exp;

	// input n: sample size
	unsigned long int n;
	cout << "Enter n: ";
	cin >> n;

	// input d: number of bins
	unsigned long int d;
	cout << "Enter d: ";
	cin >> d;

	// open file
	char filenameIn[] = "data.txt";
	ifstream fileIn(filenameIn);
	if (!fileIn.is_open())
	{
		cout << "Error opening file " << filenameIn << endl;
		return 1;
	}

	// convert sequence to [0, 1) and store in vector X
	mpf_t T[n];					 // temp vector to store mpf_t
	unsigned long int X[d] = {}; // vector to store observed frequencies
	string line;
	for (i = 0; i < n && getline(fileIn, line); i++)
	{
		mpf_init_set_str(T[i], line.c_str(), 10);
		mpf_div_2exp(T[i], T[i], m2exp);
		mpf_mul_ui(T[i], T[i], d);
		mpf_floor(T[i], T[i]);
		X[mpf_get_ui(T[i])]++;
		mpf_clear(T[i]);
	}

	// close file
	if (fileIn.bad())
	{
		perror("error while reading file");
		return 1;
	}
	fileIn.close();

	// calculating V: the ChiSquare statistic
	double V = 0;
	for (i = 0; i < d; i++)
	{
		cout << i << "    " << X[i] << endl;
		V += (double)d / n * pow(X[i] - (double)n / d, 2);
	}
	cout << "V = " << V << endl;

	// calculating the p-value
	double p = boost::math::gamma_p((d - 1) / 2, V);
	cout << "p-value = " << p * 100 << '%' << endl;

	cout << "---the end---" << endl;
	return 0;
}
