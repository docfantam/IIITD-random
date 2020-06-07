#include <iostream>
#include <stdio.h>
#include <fstream>
#include <gmp.h>
#include <cmath>

using namespace std;

/* Evaluating Kolmogorov's distribution */
void mMultiply(double *A, double *B, double *C, int m)
{
	int i, j, k;
	double s;
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
		{
			s = 0.;
			for (k = 0; k < m; k++)
				s += A[i * m + k] * B[k * m + j];
			C[i * m + j] = s;
		}
}

void mPower(double *A, int eA, double *V, int *eV, int m, int n)
{
	double *B;
	int eB, i;
	if (n == 1)
	{
		for (i = 0; i < m * m; i++)
			V[i] = A[i];
		*eV = eA;
		return;
	}
	mPower(A, eA, V, eV, m, n / 2);
	B = (double *)malloc((m * m) * sizeof(double));
	mMultiply(V, V, B, m);
	eB = 2 * (*eV);
	if (n % 2 == 0)
	{
		for (i = 0; i < m * m; i++)
			V[i] = B[i];
		*eV = eB;
	}
	else
	{
		mMultiply(A, B, V, m);
		*eV = eA + eB;
	}
	if (V[(m / 2) * m + (m / 2)] > 1e140)
	{
		for (i = 0; i < m * m; i++)
			V[i] = V[i] * 1e-140;
		*eV += 140;
	}
	free(B);
}

double K(int n, double d)
{
	int k, m, i, j, g, eH, eQ;
	double h, s, *H, *Q;
	//OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL
	s = d * d * n;
	if (s > 7.24 || (s > 3.76 && n > 99))
		return 1 - 2 * exp(-(2.000071 + .331 / sqrt(n) + 1.409 / n) * s);
	k = (int)(n * d) + 1;
	m = 2 * k - 1;
	h = k - n * d;
	H = (double *)malloc((m * m) * sizeof(double));
	Q = (double *)malloc((m * m) * sizeof(double));
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
			if (i - j + 1 < 0)
				H[i * m + j] = 0;
			else
				H[i * m + j] = 1;
	for (i = 0; i < m; i++)
	{
		H[i * m] -= pow(h, i + 1);
		H[(m - 1) * m + i] -= pow(h, (m - i));
	}
	H[(m - 1) * m] += (2 * h - 1 > 0 ? pow(2 * h - 1, m) : 0);
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
			if (i - j + 1 > 0)
				for (g = 1; g <= i - j + 1; g++)
					H[i * m + j] /= g;
	eH = 0;
	mPower(H, eH, Q, &eQ, m, n);
	s = Q[(k - 1) * m + k - 1];
	for (i = 1; i <= n; i++)
	{
		s = s * i / n;
		if (s < 1e-140)
		{
			s *= 1e140;
			eQ -= 140;
		}
	}
	s *= pow(10., eQ);
	free(H);
	free(Q);
	return s;
}

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

	// open file
	char filenameIn[] = "data.txt";
	ifstream fileIn(filenameIn);
	if (!fileIn.is_open())
	{
		cout << "Error opening file " << filenameIn << endl;
		return 1;
	}

	// convert sequence to [0, 1) and store in vector X
	mpf_t T[n];  // temp vector to store mpf_t
	double X[n]; // vector to store [double] inputs
	string line;
	for (i = 0; i < n && getline(fileIn, line); i++)
	{
		mpf_init_set_str(T[i], line.c_str(), 10);
		mpf_div_2exp(T[i], T[i], m2exp);
		X[i] = mpf_get_d(T[i]); // X[i] belongs to [0, 1)
		cout << X[i] << endl;
		mpf_clear(T[i]);
	}

	// close file
	if (fileIn.bad())
	{
		perror("error while reading file");
		return 1;
	}
	fileIn.close();

	// step 1 [Initialize]
	double MIN[n + 1];
	for (i = 0; i <= n; i++)
	{
		MIN[i] = 1;
	}
	double MAX[n + 1] = {};
	double NUM[n + 1] = {};

	// step 2 [Put input observations into bins]
	double f;
	unsigned long int j;

	for (i = 1; i <= n; i++)
	{
		f = X[i - 1];
		j = ceil(f * n); // compute bin for X
		NUM[j]++;
		if (MAX[j] < f)
			MAX[j] = f;
		if (MIN[j] > f)
			MIN[j] = f;
	}

	// step 3 [Process each bin; find maximum positive and negative deviates
	j = 0;
	double DP = 0;
	double DN = 0;
	double z;

	for (i = 0; i <= n; i++)
	{
		if (NUM[i] > 0)
		{
			z = MIN[i] - (double)j / n;
			if (z > DN)
				DN = z;
			j = j + NUM[i];
			z = (double)j / n - MAX[i];
			if (z > DP)
				DP = z;
		}
	}

	// step 4 [Compute KPmax, KNmax, and Kmax]
	double KPmax = sqrt(n) * DP;
	double KNmax = sqrt(n) * DN;
	double Kmax = max(KPmax, KNmax);

	cout << "K+max = " << KPmax << endl;
	cout << "K-max = " << KNmax << endl;
	cout << "Kmax = " << Kmax << endl;

	cout << "p-value = " << K(n, Kmax) * 100 << '%' << endl;
	cout << "---the end---" << endl;
	return 0;
}
