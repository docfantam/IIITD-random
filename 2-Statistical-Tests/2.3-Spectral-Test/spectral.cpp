#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <getopt.h>

#include "spectral.hpp"

using namespace std;

int main(int argc, char *argv[])
{
	int option;
	bool verbose = false;
	while ((option = getopt(argc, argv, "v")) != -1)
	{
		switch (option)
		{
		case 'v':
			verbose = 1;
			break;
		default:
			std::cerr << "Usage: " << argv[0] << " [-v Verbose]" << std::endl;
			return 1;
			break;
		}
	}
	spectral(verbose);
	std::cout << "---the end---" << std::endl;
	return 0;
}
