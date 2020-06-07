#include <iostream>
#include <stdio.h>
#include <gmp.h>

#define T 9

using namespace std;

void disp(mpz_t U[T + 1][T + 1], int t)
{
	for (int i = 1; i <= t; i++)
	{
		for (int j = 1; j <= t; j++)
		{
			gmp_printf("%12Zd\t", U[i][j]);
		}
		std::cout << std::endl;
	}
}

void dotprod(mpz_t product, mpz_t U[T + 1], mpz_t V[T + 1], int t)
{
	mpz_set_ui(product, 0);
	for (int i = 1; i <= t; i++)
	{
		mpz_addmul(product, U[i], V[i]);
	}
}

void divround(mpz_t q, mpz_t n, mpz_t d)
{
	mpf_t temp1;
	mpf_init(temp1);
	mpf_set_z(temp1, n);
	mpf_t temp2;
	mpf_init(temp2);
	mpf_set_z(temp2, d);
	mpf_div(temp1, temp1, temp2);
	mpf_set_d(temp2, 0.5);
	mpf_add(temp1, temp1, temp2);
	mpf_floor(temp1, temp1);
	mpz_set_f(q, temp1);
}

void spectral(bool verbose)
{
	char inputStr[1024];

	mpz_t temp1, temp2, temp3; // temp variables
	int i, j, k, l;			   // iterators
	int count;
	mpz_init(temp1);
	mpz_init(temp2);
	mpz_init(temp3);

	mpz_t X[T + 1], Y[T + 1], Z[T + 1];
	for (i = 1; i <= T; i++)
	{
		mpz_init(X[i]);
		mpz_init(Y[i]);
		mpz_init(Z[i]);
	}

	// step 0
	mpz_t a, c, m;
	std::cout << "Enter a: ";
	std::cin >> inputStr;
	mpz_init_set_str(a, inputStr, 10);
	std::cout << "Enter c: ";
	std::cin >> inputStr;
	mpz_init_set_str(c, inputStr, 10);
	std::cout << "Enter m: ";
	std::cin >> inputStr;
	mpz_init_set_str(m, inputStr, 10);

	// step 1
	int t = 2;
	mpz_t h, h_, p, p_, r, s;
	mpz_init_set(h, a);		// h = a
	mpz_init_set(h_, m);	// h_ = m
	mpz_init_set_si(p, 1);	// p = 1
	mpz_init_set_si(p_, 0); // p = 0
	mpz_init_set(r, a);		// r = a
	mpz_init_set_ui(s, 1);
	mpz_addmul(s, a, a); // s = 1 + a*a

	// step2
S2:
	mpz_t q, u, v;

	mpz_init(q);
	mpz_fdiv_q(q, h_, h); // q = h_ / h

	mpz_init_set(u, h_);
	mpz_submul(u, q, h); // u = h_ - q*h

	mpz_init_set(v, p_);
	mpz_submul(v, q, p); // v = p_ - q*p

	mpz_mul(temp1, u, u);
	mpz_mul(temp2, v, v);
	mpz_add(temp3, temp1, temp2);
	if (mpz_cmp(s, temp3) > 0)
	{					   // if u*u + v*v < s
		mpz_set(s, temp3); // s = u*u + v*v
		mpz_set(h_, h);	   // h_ = h
		mpz_set(h, u);	   // h = u
		mpz_set(p_, p);	   // p_ = p
		mpz_set(p, v);	   // p = v
		goto S2;		   // repeat S2
	}

	// step 3
	mpz_sub(u, u, h); // u = u - h
	mpz_sub(v, v, p); // v = v - p
	mpz_mul(temp1, u, u);
	mpz_mul(temp2, v, v);
	mpz_add(temp3, temp1, temp2);
	if (mpz_cmp(s, temp3) > 0)
	{					   // if u*u + v*v < s
		mpz_set(s, temp3); // s = u*u + v*v
		mpz_set(h_, u);	   // h_ = u
		mpz_set(p_, v);	   // p_ = v
	}
	mpz_sqrt(temp1, s);
	gmp_printf("v2 = %Zd\n", temp1); // echo v2 = sqrt(s)

	mpz_t U[T + 1][T + 1], V[T + 1][T + 1];
	for (i = 0; i <= T; i++)
	{
		for (j = 0; j <= T; j++)
		{
			mpz_init(U[i][j]);
			mpz_init(V[i][j]);
		}
	}

	mpz_neg(temp1, h);
	mpz_neg(temp2, h_);
	mpz_neg(temp3, p);

	mpz_set(U[1][1], temp1); // U[1][1] = -h
	mpz_set(U[1][2], p);	 // U[1][2] = p
	mpz_set(U[2][1], temp2); // U[2][1] = -h_
	mpz_set(U[2][2], p_);	 // U[2][2] = p_

	mpz_set(V[1][1], p_);	 // V[1][1] = p_
	mpz_set(V[1][2], h_);	 // V[1][2] = h_
	mpz_set(V[2][1], temp3); // V[2][1] = -p
	mpz_set(V[2][2], temp1); // V[2][2] = -h

	if (mpz_sgn(p_) == 1)
	{ // if p_ > 0
		for (i = 1; i <= 2; i++)
		{ // V = -V
			for (j = 1; j <= 2; j++)
			{
				mpz_neg(V[i][j], V[i][j]);
			}
		}
	}

	// step 4
S4:
	if (t == T)
	{
		return;
	}
	t++; // t = t+1
	mpz_mul(r, a, r);
	mpz_fdiv_r(r, r, m); // r = (a*r) mod m
	mpz_neg(temp1, r);
	mpz_set(U[t][1], temp1); // U[t][1] = -r
	mpz_set_si(U[t][t], 1);	 // U[t][t] = 1
	mpz_set(V[t][t], m);	 // V[t][t] = m

	for (i = 1; i < t; i++)
	{ // for 1 <= i < t
		mpz_mul(temp1, V[i][1], r);
		//mpz_fdiv_q(q, temp1, m);		// q = V[i][1]*r / m
		divround(q, temp1, m);
		mpz_mul(temp2, q, m);
		mpz_sub(V[i][t], temp1, temp2); // V[i][t] = V[i][1]*r - q*m
		for (j = 1; j <= t; j++)
		{
			mpz_addmul(U[t][j], q, U[i][j]); // U[t] = U[t] + q*U[i]
		}
	}

	dotprod(temp1, U[t], U[t], t);
	if (mpz_cmp(s, temp1) > 0)
	{
		mpz_set(s, temp1); // s = min(s, U[t].U[t])
	}

	k = t; // k = t
	j = 1; // j = 1

	if (verbose)
	{
		std::cout << "-----After step 4-----" << std::endl;
		std::cout << "U: \n";
		disp(U, t);
		std::cout << "V: \n";
		disp(V, t);
	}
	count = 0;

	// step 5
S5:
	if (verbose)
		std::cout << "-----Entered step 5-----" << std::endl;

	count++;
	for (i = 1; i <= t; i++)
	{
		dotprod(temp1, V[i], V[j], t); // temp1 = V[i].V[j]
		mpz_abs(temp3, temp1);
		mpz_mul_si(temp3, temp1, 2);				 // temp3 = 2|V[i].V[j]|
		dotprod(temp2, V[j], V[j], t);				 // temp2 = V[j].V[j]
		if ((i != j) && (mpz_cmp(temp3, temp2) > 0)) // if 2|V[i].V[j]| > V[j].V[j]
		{
			divround(q, temp1, temp2); // q = V[i].V[j] / V[j].V[j]
			for (l = 1; l <= t; l++)
			{
				mpz_submul(V[i][l], q, V[j][l]); // V[i] = V[i] - q*V[j]
				mpz_addmul(U[j][l], q, U[i][l]); // U[j] = U[j] + q*U[i]
			}
			dotprod(temp3, U[j], U[j], t);
			if (mpz_cmp(s, temp3) > 0)
			{
				mpz_set(s, temp3); // s = min(s, U[j].U[j]
			}
			k = j; // k = j
		}
		if (verbose)
		{
			std::cout << "count: " << count << "    ";
			std::cout << "i: " << i << "     ";
			gmp_printf("q: %Zd\n", q);
			std::cout << "U: \n";
			disp(U, t);
			std::cout << "V: \n";
			disp(V, t);
			std::cout << "--------------------------------------------------------------------------------" << std::endl;
		}
	}

	// step 6
	j = (j == t) ? 1 : j + 1; // if j=t, set j=1; otherwise set j=j+1
	if (j != k)				  // if j!=k, return to S5
	{
		goto S5;
	}

	// step 7
	if (verbose)
		std::cout << "-----Entered step 7-----" << std::endl;

	for (l = 1; l <= t; l++)
	{
		mpz_set_si(X[l], 0); // X <- (0, ..., 0)
		mpz_set_si(Y[l], 0); // Y <- (0, ..., 0)
	}
	k = t;					 // k = t
	for (j = 1; j <= t; j++) // for 1 <= j <= t
	{
		dotprod(temp1, V[j], V[j], t);
		mpz_mul(temp1, temp1, s);
		mpz_mul(temp2, m, m);
		mpz_fdiv_q(temp1, temp1, temp2);
		// gmp_printf("temp1 = %Zd\ntemp2 = %Zd\n", temp1, temp2);
		mpz_sqrt(Z[j], temp1); // Z[j] = sqrt(V[j].V[j] * s / m / m)
	}

	// step 8
S8:
	if (mpz_cmp(X[k], Z[k]) == 0)
	{
		goto S10;
	}
	else
	{
		mpz_add_ui(X[k], X[k], 1);
		for (l = 1; l <= t; l++)
		{
			mpz_add(Y[l], Y[l], U[k][l]);
		}
	}

	// step 9
S9:
	k++;
	if (k <= t)
	{						 // if k <= t
		mpz_neg(X[k], Z[k]); // X[k] = -Z[k]
		for (l = 1; l <= t; l++)
		{
			mpz_mul(temp1, Z[k], U[k][l]);
			mpz_submul_ui(Y[l], temp1, 2); // Y = Y - 2*Z[k]*U[k]
		}
		goto S9;
	}
	else
	{ // but if k > t
		dotprod(temp1, Y, Y, t);
		if (mpz_cmp(s, temp1) > 0)
		{ // s = min(s, Y.Y)
			mpz_set(s, temp1);
		}
	}

	// step 10
S10:
	k--;
	if (k >= 1)
	{
		goto S8;
	}
	else
	{
		mpz_sqrt(temp1, s);
		gmp_printf("v%d = %Zd\n", t, temp1); // echo v_t = sqrt(s)
		goto S4;
	}
}
