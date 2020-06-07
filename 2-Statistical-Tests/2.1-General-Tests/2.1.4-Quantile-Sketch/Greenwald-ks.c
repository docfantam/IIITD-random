//*************** Include files ***********************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gmp.h>

//*************** Structure Declarations **************************************
struct Summary {
	mpz_t value;
	long int rank;
	long int wiggle;
	struct Summary * next;
	struct Summary * prev;
};

//*************** Constant Definitions ****************************************
#define TRUE 1
#define FALSE 0
typedef struct Summary Summary;

//*************** Function Declarations ***************************************
long int band();
void Insert();
void COMPRESS();
void PrintCount();
void PrintStructure();
void ks();
//*************** Global Variable Declarations ********************************
long int window;
long int counter;
int useband;
int verbose;
Summary * S = NULL;
double epsilon;
long int datapoints = 100000;
Summary * current = NULL;

//*************** Here starts the main program ********************************

int main ( int argc, char * argv[] )
{ //--> BEGIN MAIN PROGRAM

	// Variable Declaration & Initialization
	// INPUT VARIABLES
	mpz_t dp ;
	mpz_init(dp);

	// OUTPUT VARIABLES

	// PROCESSING VARIABLES
	long int i ;
	long int delta ;


	/* HERE START THE EXECUTABLE STATEMENTS */

	// get command line arguments
	epsilon = atof ( argv[2] ) ;
	useband = ( argv[3][0] == 'u' ) ;
	verbose = ( argv[4][0] == 'v' ) ;	

	/* MT */
	char inputStr[1024];

    // Random State Initialization [see 9.1]
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    // Random State Seeding [see 9.2]
    mpz_t seed;
    printf("Enter seed: ");
    scanf("%s", inputStr);
    mpz_init_set_str(seed, inputStr, 10);
    gmp_randseed(state, seed);
    
    mp_bitcnt_t m2exp = 64;
 
	counter = 0 ;
	window = 1 / ( 2 * epsilon ) ;
	for ( i = 0; i < datapoints; i++ )
	{
		if ( ( counter % window ) == 0 )
		{
			if ( verbose )
			{
				PrintStructure();
				COMPRESS();
				PrintStructure();
			}
			else
			{
				COMPRESS();
			}			
		}

		mpz_urandomb(dp, state, m2exp);

		Insert(dp);

		counter++ ;
	}

	printf ( "FINAL RESULT:\n" ) ;

	PrintCount();

	PrintStructure();

	ks();

	// end the program
	return 0;

} //--> END MAIN PROGRAM


//*************** Here are the function definitions ***************************
void PrintStructure()
{
	Summary * out ;

	if ( S == NULL ) return;

	out = S ;

	printf ( "S: " ) ;

	while ( out != NULL )
	{
		gmp_printf ( "(%Zd,%ld,%ld) ", out->value, out->rank, out->wiggle ) ;

		out = out->next ;
	}

	printf ( "\n\n" ) ;
}

void PrintCount()
{
	long int count = 0 ;
	Summary * temp ;

	temp = S ;

	while ( temp != NULL )
	{
		count++;
		temp = temp->next ;
	}

	printf ( "Count: %d\n", count ) ;
}

long int FindEstimate(mpz_t value)
{
	long int rank = 0;
	long int margin ;
	Summary * temp ;

	temp = S;

	while ( temp != NULL )
	{
		if ( mpz_cmp(temp->value, value)>0 )
		{
			return rank;
		}
		rank += temp->rank ;
		temp = temp->next ;
	}
}

void COMPRESS()
{
	long int valuecomp ;
	long int comparevalue ;
	long int bandone, bandtwo ;

	if ( (S == NULL) || (S->next == NULL) )
		return;

	comparevalue = 2 * epsilon * counter ;

	current = S->next ;

	while ( current != NULL )
	{
		if (current->next == NULL )
			break ;

		current = current->next ;
	}

	current = current->prev ;

	while ( current != NULL )
	{
		bandone = band ( current->wiggle, comparevalue ) ;
		bandtwo = band ( current->next->wiggle, comparevalue ) ;
		valuecomp = current->rank + current->next->rank + current->next->wiggle ;

		if ( (bandone <= bandtwo) && (valuecomp < comparevalue) )
		{
			current->next->rank += current->rank ;

			if ( current->prev != NULL )
			{
				current->prev->next = current->next ;
			}
			else
			{
				S = current->next ;
			}

			current->next->prev = current->prev ;
		}
		current = current->prev ;
	}
}

long int band (long int p, long int second )
{
	long int band = 0;
	long int pow2 = 1;
	long int i, k1, k2;

	if ( ! useband )
	{
		return 0 ;
	}

	k1 = p - pow2 - p%pow2;

	if ( p == second )
		return 0 ;

	while ( k1 >= 0 )
	{
		k2 = k1;

		band++;

		pow2 = 2*pow2;

		k1 = p - pow2 - p%pow2;

		if ( k1 >= 0 )
		{
			for (i = k1+1; i <= k2; i++)
				if ( second == i )
					return band ;
		}
	}

	if ( second == 0 )
		return band ;

	return 0 ;
}

void Insert ( mpz_t val )
{	
	Summary * new = (Summary*) malloc ( sizeof(Summary) ) ;
	mpz_init_set(new->value, val);
	new->rank = 1;
	new->next = NULL;
	new->prev = NULL;
	
	// empty
	if(S == NULL)
	{
		new->wiggle = 0;
		S = new;
		return;
	}
	
	// non empty
	Summary * previous = NULL;
	current = S ;
	while ( current != NULL )
	{
		if ( mpz_cmp(current->value, val) > 0 )
		{
			break ;
		}

		previous = current ;
		current = current->next ;
	}
	// if insertion at beginning
	if(previous == NULL)
	{
		new->wiggle = 0;
		new->next = current;
		current->prev = new;
		S = new;
	}
	// if insertion at end
	else if(current == NULL)
	{
		new->wiggle = 0;
		previous->next = new;
		new->prev = previous;
	}
	// if somewhere in the middle
	else
	{
		new->wiggle = (int) ( 2 * epsilon * counter );
		new->prev = previous;
		previous->next = new;
		new->next = current;
		current->prev = new;
	}
}

void ks()
{
	/* Find size of Summary */
	long int count = 0 ;
	Summary * temp = S ;
	mpz_t xZ; mpz_init(xZ);
	mpf_t xF; mpf_init(xF);
	double E;
	double D = 0;
	
	while ( temp->next != NULL ) // PAY ATTENTION: dont go till temp != NULL because the last node has the remaining rank (large)
	{
		count++;
		
		mpz_set(xZ, temp->value);
		mpf_set_z(xF, xZ);
		mpf_div_2exp(xF, xF, 64);
		//gmp_printf("%Ff\t", xF);
		E = fabs((double)FindEstimate(xZ)/(double)datapoints - mpf_get_d(xF));
		printf("%f\t", E);
		if(D < E)
			D = E;
		
		temp = temp->next ;
	}
	printf("KS statistic = %lf\n", D);
}
	
		
