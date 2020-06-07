#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <cmath>                                  // for exp(), sqrt(), pow()
#include <boost/math/special_functions/gamma.hpp> // for normalized incomplete gamma functionI

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
start:
    // Choosing RNG
    char rng;
    cout << "Enter choice of RNG: " << endl;
    cout << "\t(a) LCG (Linear Congruential Generator)" << endl;
    cout << "\t(b) MT19937 (Mersenne Twister)" << endl;
    cin >> rng;

    // Random State
    gmp_randstate_t state;

    char inputStr[1024];
    mpz_t seed;

    // LCG inputs
    mpz_t a;           // a: multiplier
    unsigned long c;   // c: increment
    mp_bitcnt_t m2exp; // 2^m2exp: modulus

    switch (rng)
    {
    case 'a': // LCG
        // input a
        cout << "Enter a: ";
        cin >> inputStr;
        mpz_init_set_str(a, inputStr, 10);

        // input c
        cout << "Enter c: ";
        cin >> c;

        // input m2exp
        cout << "Enter m2exp: ";
        cin >> m2exp;

        // Random State Initialization [see 9.1]
        gmp_randinit_lc_2exp(state, a, c, m2exp);

        // Random State Seeding [see 9.2]
        cout << "Enter seed: ";
        cin >> inputStr;
        mpz_init_set_str(seed, inputStr, 10);
        gmp_randseed(state, seed);
        break;

    case 'b': //MT
        // Random State Initialization [see 9.1]
        gmp_randinit_mt(state);

        // Random State Seeding [see 9.2]
        cout << "Enter seed: ";
        cin >> inputStr;
        mpz_init_set_str(seed, inputStr, 10);
        gmp_randseed(state, seed);
        break;
    default:
        goto start;
    }

    //---------------------------------------------
    unsigned long int i; // loop iterator

    // input n: sample size (common for both tests)
    unsigned long int n;
    cout << "Enter n: ";
    cin >> n;

    // input d: number of bins (for ChiSquare test)
    unsigned long int d;
    cout << "Enter d: ";
    cin >> d;

    // step 1 [Initialize]
    double MIN[n + 1];
    for (i = 0; i <= n; i++)
    {
        MIN[i] = 1;
    }
    double MAX[n + 1] = {};
    double NUM[n + 1] = {};

    // step 2 [Put input observations into bins]

    /* we need to convert the numbers into [0, 1)
       in case of LCG, we divide by 2^m2exp
       in case of MT,  we divide by 2^word_length
    */
    if (rng == 'b') // MT
    {
        m2exp = 64; // word length
    }

    mpz_t xZ; // holds the generated RNs (integer)
    mpz_init(xZ);

    mpf_t xF; // holds the generated RNs (converted to floating point)
    mpf_init(xF);

    unsigned long int X[d] = {}; // vector to store observed frequencies for ChiSquare test

    double f;
    unsigned long int j;

    for (i = 1; i <= n; i++)
    {
        mpz_urandomb(xZ, state, m2exp); // generate next RN
        mpf_set_z(xF, xZ);
        mpf_div_2exp(xF, xF, m2exp); // normalize to [0, 1)
        f = mpf_get_d(xF);           // convert to double [Note: F(X) = X here]
        j = ceil(f * n);             // compute bin for KS test
        mpf_mul_ui(xF, xF, d);
        mpf_floor(xF, xF);   // compute bin for ChiSquare test
        X[mpf_get_ui(xF)]++; // update frequency
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
    double p; // p-value

    // calculating V: the ChiSquare statistic
    double V = 0;
    for (i = 0; i < d; i++)
    {
        V += (double)d / n * pow(X[i] - (double)n / d, 2);
    }

    // Output results for KS test
    cout << "Results for KS test: " << endl;
    cout << "\tK+max = " << KPmax << endl;
    cout << "\tK-max = " << KNmax << endl;
    cout << "\tKmax = " << Kmax << endl;
    p = K(n, Kmax) * 100;
    cout << "\tp-value = " << p << '%' << endl
         << endl;

    // Output results for ChiSquare test
    cout << "Results for ChiSquare test: " << endl;
    cout << "\tV = " << V << endl;
    p = boost::math::gamma_p((d - 1) / 2, V) * 100;
    cout << "\tp-value = " << p << '%' << endl;

    cout << "---the end---" << endl;
    return 0;
}
