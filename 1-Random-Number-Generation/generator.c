#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <unistd.h>

#include "LCG.h"
#include "MT.h"

int main(int argc, char *argv[])
{
	/* ---for getopt--- */
	extern char* optarg;
	extern int optind;

	char usage[100] = "usage: %s -h | -g gname -p pname [-n nname] [-f fname] \n";

	/* ---option read--- */
	int c;

	/* ---option flags--- */
	int hflag = 0;	// help
	int gflag = 0;	// generator
	int pflag = 0;	// parameter
	int nflag = 0;	// size
	int fflag = 0;	// file
	int eflag = 0;	// error

	/* ---option names--- */
	char *gname = NULL;	// generator
	char *pname = NULL;	// parameter
	char *nname = NULL;	// size
	char *fname = NULL;	// file

	while((c = getopt(argc, argv, "hg:p:n:f:e")) != -1)
	{
		switch(c)
		{
			case 'h': hflag = 1; break;
			case 'g': gflag = 1; gname = optarg; printf("gname:%s\n", gname); break;
			case 'p': pflag = 1; pname = optarg; printf("pname:%s\n", pname); break;
			case 'n': nflag = 1; nname = optarg; printf("nname:%s\n", nname); break;
			case 'f': fflag = 1; fname = optarg; printf("fname:%s\n", fname); break;
			case '?': eflag = 1; break;
		}
	}

	if(hflag == 1) 								// help
	{
		printf("%s\n", __rng_desc_lcg());
		printf("%s\n", __rng_desc_mt());
	}
	else if(gflag == 1)							// generator
	{
		if(strcmp(gname, __rng_name_lcg()) == 0)
		{
			__rng_run_lcg(pname, fname, nname);
		}
		else if(strcmp(gname, __rng_name_mt()) == 0)
		{
			__rng_run_mt(pname, fname, nname);
		}	
	}
	else if(eflag == 1)							// error
	{
		fprintf(stderr, usage, argv[0]);
		exit(1);
	}

	printf("---the end---\n");
	return 0;
}