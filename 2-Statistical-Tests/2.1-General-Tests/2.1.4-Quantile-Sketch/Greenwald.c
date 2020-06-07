/* From rlohfer: file `/home/rlohfer/cs584/Greenwald.c' as assignment hw5-0 */
/*
THIS CODE IS MY OWN WORK, IT WAS WRITTEN WITHOUT CONSULTING 
A TUTOR OR CODE WRITTEN BY OTHER STUDENTS -
Robin Lohfert
*/

/* File Header:
ASSIGNMENT:	
DIRECTORY:	Network
FILENAME:	
PROGRAMMER:	Robin Lohfert
IDE:		XCode & geedit 
COMPILER:	GCC & CC
RELATED FILES:	
DATE:		4/14/2006
DATE DUE:	4/20/2006
INSTRUCTOR:	Dr. Cheung
CLASS:		CS 584
*/

//*************** Include files ***********************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//*************** Structure Declarations **************************************
struct Summary {
	int value;
	int rank;
	int wiggle;
	struct Summary * next;
	struct Summary * prev;
};

//*************** Constant Definitions ****************************************
#define TRUE 1
#define FALSE 0
typedef struct Summary Summary;

//*************** Function Declarations ***************************************
int band();
void Insert();
void COMPRESS();
void PrintCount();
void FindPosition();
void PrintResults();
void PrintStructure();

//*************** Global Variable Declarations ********************************
FILE *data;
int window;
int * Input;
int counter;
int useband;
int verbose;
Summary * S;
double epsilon;
int datapoints;
int HeadOrTail;
Summary * current;

int intcompare (const int *a, const int *b)
{
	return (int) ( *a - *b ) ;
}
//*************** Here starts the main program ********************************

int main ( int argc, char * argv[] )
{ //--> BEGIN MAIN PROGRAM

	// Variable Declaration & Initialization
	// INPUT VARIABLES
	int dp ;

	// OUTPUT VARIABLES

	// PROCESSING VARIABLES
	int i ;
	int delta ;


	/* HERE START THE EXECUTABLE STATEMENTS */


	// open the input data file
	data = fopen ( argv[1], "r" ) ;

	// get command line arguments
	epsilon = atof ( argv[2] ) ;
	useband = ( argv[3][0] == 'u' ) ;
	verbose = ( argv[4][0] == 'v' ) ;

	// Read N data, get rid of it
	fscanf ( data, "%d", &datapoints ) ;

	// allocate space to remember input
	Input = (int*) malloc ( datapoints * sizeof(int) ) ; 


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

		fscanf ( data, "%d", &dp ) ;

		Input[i] = dp ;

		FindPosition ( dp ) ;

		if ( HeadOrTail )
		{
			delta = 0 ;
		}
		else
		{
			delta = (int) ( 2 * epsilon * counter ) ;
		}

		Insert ( dp, 1, delta ) ;

		counter++ ;
	}

	qsort ( Input, datapoints, sizeof(int), (void*)intcompare ) ;

	printf ( "FINAL RESULT:\n" ) ;

	PrintCount();

	PrintStructure();

	PrintResults();


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
		printf ( "(%d,%d,%d) ", out->value, out->rank, out->wiggle ) ;

		out = out->next ;
	}

	printf ( "\n\n" ) ;
}

void PrintCount()
{
	int count = 0 ;
	Summary * temp ;

	temp = S ;

	while ( temp != NULL )
	{
		count++;
		temp = temp->next ;
	}

	printf ( "Count: %d\n", count ) ;
}

int FindEstimate(int i)
{
	int rank ;
	int margin ;
	Summary * temp ;

	temp = S ;
	rank = 0 ;
	margin = epsilon * counter ;

	while ( temp != NULL )
	{
		rank += temp->rank ;

		if ( rank >= i+1 )
			return temp->value;

		temp = temp->next ;
	}

	return -1 ;
}

int FindError(int start, int actual, int estimate )
{
	int i, margin ;

	margin = epsilon * counter ;

	for ( i = 0; i < datapoints; i++ )
	{
		if ( start+i < datapoints )
			if ( Input[start+i] == estimate )
				return i ;

		if ( start-i >= 0 )
			if ( Input[start-i] == estimate )
				return i ;
	}

	return 999 ;
}

void PrintResults()
{
	int estimate, error, i ;

	printf ( "%.2f\n", epsilon * counter ) ;

	for ( i = 0; i < datapoints; i++ )
	{
		estimate = FindEstimate(i);
		error = FindError ( i, Input[i], estimate ) ;

		printf ( "%d\t", i+1 ) ;
		printf ( "%d\t", Input[i] ) ;
		printf ( "%d\t", estimate ) ;
		printf ( "%d\n", error ) ;
	}
}

void COMPRESS()
{
	int valuecomp ;
	int comparevalue ;
	int bandone, bandtwo ;

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

int band ( int p, int second )
{
	int band = 0;
	int pow2 = 1;
	int i, k1, k2;

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

void Insert ( int value, int rank, int delta )
{
	Summary * new ;

	new = (Summary*) malloc ( sizeof(Summary) ) ;

	new->value = value ;
	new->rank = rank ;
	new->wiggle = delta ;

	if ( current == NULL )
	{
		new->next = S ;
		new->prev = NULL ;
		S = new ;
	}
	else
	{
		new->next = current->next ;
		new->prev = current ;
		current->next = new ;
	}

	if ( new->next != NULL )
		new->next->prev = new ;
}

void FindPosition ( int val )
{
	Summary * previous ;

	current = S ;
	previous = NULL ;

	while ( current != NULL )
	{
		if ( current->value > val )
		{
			break ;
		}

		previous = current ;
		current = current->next ;
	}

	if ( previous == NULL )
	{
		HeadOrTail = TRUE ;
	}
	else if ( current == S )
	{
		HeadOrTail = TRUE ;
		current = NULL ;
	}
	else if ( current == NULL )
	{
		current = previous ;
		HeadOrTail = TRUE ;
	}
	else
	{
		current = previous ;
		HeadOrTail = FALSE ;
	}
}
