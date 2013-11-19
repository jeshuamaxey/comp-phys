#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include "dimensions.cpp"

#define simulatedTime 100.0												//simulated time (seconds)
#define g 9.81																		//acceleration due to gravity
#define pi atan(1.0)															//mutha fuckin pi man

#define N 10																			//there are NxN spins simulated

/****** PROGRAM CONTROL SETTINGS ******/
#define coldStart false														//use for starting at low T

using namespace std;

/****** FUNCTION PROTOTYPES ***********/
gsl_rng* setupUniformRNG();
void initialiseSpins();
double calcTotalEnergy();
void findEquilibrium(double);
void outputSpinsToFile();
void alignSpins();
void randomlyDistSpins();
int randomSpin();
int randIndex(int);

void updateProgress(double);
void done();

int spin[N][N];
double total_E;
double J = 1.0;
//for gsl rng help see here:
//http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-initialization.html#Random-number-generator-initialization
gsl_rng * r_uni;						//default rng instance

int main()
{
	double initial_temp = 0.0;

	
	r_uni = setupUniformRNG();
	initialiseSpins();
	total_E = calcTotalEnergy();
	findEquilibrium(initial_temp);

	/*********************************/
	long int longtime = 10;
	for (int c = 0; c < longtime; ++c)
	{
		//updateProgress((1.0*c/longtime));
	}
	/*********************************/
	outputSpinsToFile();
	done();
}

/* START ALL FUNCTIONS CALLED FROM main() */

gsl_rng* setupUniformRNG()
{
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	unsigned long seed = 116426264;
	gsl_rng_set(r,seed);
	return r;
}

void initialiseSpins()
{
	if(coldStart)
	{
		alignSpins();
	}
	else
	{
		randomlyDistSpins();
	}
}

double calcTotalEnergy()
{
	double E = 0.0;
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			//the modulo division ensures the periodic boundary conditions are met
			// E += spin[x,y]*spin[(x+1.0)%N, y];
			// E += spin[x,y]*spin[x, (y+1.0)%N];
			if(x!=N-2)	E += spin[x,y]*spin[x+1,y];
			else 				E += spin[x,y]*spin[0, y];
			if(y!=N-2)	E += spin[x,y]*spin[x,y+1];
			else 				E += spin[x,y]*spin[x, 0];
		}
	}
	E *= -0.5*J;
	cout << E << "\n";
	return E;
}

void findEquilibrium(double initial_temp)
{
	bool equilibrium = false;						//used to keep track of state of equilibrium
	int s;															//the spin value for a particular loop iteration
	double dE;													//the change in energy caused by flipping the spin s

	double	E, M,												//total energy of system, magnetisation of system
					E_av, E_s_av,								//average energy, average square energy
					S_av, S_s,av;								//average spin, average square spin
	while(!equilibrium)
	{
		s = spin[randIndex(N)][randIndex(N)];

		equilibrium = true;
	}
}

void outputSpinsToFile()
{
	ofstream spins("data/spins.tsv");
	spins << "X\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\n";
	for (int x = 0; x < N; ++x)
	{
		spins << x << "\t";
		for (int y = 0; y < N; ++y)
		{
			spins << spin[x][y] << "\t";
		}
		spins << "\n";
	}
}

/* END ALL FUNCTIONS CALLED FROM main() */

void alignSpins()
{
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			spin[x][y] = 1;
		}
	}
}

void randomlyDistSpins()
{
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			spin[x][y] = randomSpin();
		}
	}
}

int randomSpin()
{
	return gsl_rng_uniform(r_uni) > 0.5 ? 1 : -1 ;
}

int randIndex(int max)
{
	return 1;
}

//updates the progress display on commandline
void updateProgress(double progress)
{
	int width = getTerminalWidth();
	int progressBarLength = int(progress*width);
	stringstream progressBar;

	for(int c = 0; c<progressBarLength; ++c)
	{
		progressBar << "#";
	}
	std::cout << "\r" << progressBar.str() << std::flush;
}

//displays done message to terminal
void done()
{
	std::cout << std::fixed << "\r"<< "\n100\% - all done. Simulated time elapsed: " << simulatedTime << "s\n" << std::flush;
}