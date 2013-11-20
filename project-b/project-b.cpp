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
#define pi 4.0*atan(1.0)													//mutha fuckin pi man

#define N 10																			//there are NxN spins simulated
#define beta 0																		// J/(k_b*T) -> 1/(k_b*T) = J*beta

/****** PROGRAM CONTROL SETTINGS ******/
#define coldStart false														//use for starting at low T

using namespace std;

/****** FUNCTION PROTOTYPES ***********/

void findEquilibrium(double);

//spin functions
void initialiseSpins();
void outputSpinsToFile();
void alignSpins();
void randomlyDistSpins();
int randomSpin();
void flipSpin(int, int);

//random number functions
gsl_rng* setupUniformRNG();
int randInt(int);

//energy function
double calcTotalEnergy();
double calcSiteEnergy(int, int);
double calcDeltaEnergy(int, int);

//magnetisation functions
double calcTotalMagnetisation();

//helper functions
void updateProgress(double);
void done();

int spin[N][N];
double J = 1.0;
//for gsl rng help see here:
//http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-initialization.html#Random-number-generator-initialization
gsl_rng * r_uni;						//default rng instance

int main()
{
	double initial_temp = 0.0;

	
	r_uni = setupUniformRNG();
	initialiseSpins();
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

void findEquilibrium(double initial_temp)
{
	bool equilibrium = false;							//used to keep track of state of equilibrium
	int x,y;															//the spin coordinates particular loop iteration
	double dE;														//the change in energy caused by flipping the spin s

	double	E = calcTotalEnergy(),				//total energy of system
					M = calcTotalMagnetisation(),	//total magnetisation of system
					E_av, E_s_av,									//average energy, average square energy
					S_av, S_s,av;									//average spin, average square spin
	while(!equilibrium)
	{
		x = randInt(N);
		y = randInt(N);
		dE = calcDeltaEnergy(x,y);
		if(dE < 0)
		{
			flipSpin(x,y);
			E += dE;
		}
		else
		{
			double r = gsl_rng_uniform(r_uni); 			//random number 0 < r < 1
			if(r < exp(-dE*beta) )
			{
				flipSpin(x,y);
				E += dE;
			}
		}

		equilibrium = true;
	}
}

/*
* SPIN MANIPULATION
*/

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

void flipSpin(int x, int y)
{
	spin[x][y] = -1*spin[x][y];
}

/*
*	ENERGY FUNCTIONS
*/

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

int randInt(int max)
{
	return floor( gsl_rng_uniform(r_uni)*max );
}

/*
*	ENERGY FUNCTIONS
*/
double calcTotalEnergy()
{
	double E = 0.0;
	int right, down, left, up;
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			E += 0.5*calcSiteEnergy(x,y);
		}
	}
	//multiply by factor at front of equation at end because I am an efficient coder LOLJK
	cout << "Energy: " << E << "\n";
	return E;
}

double calcSiteEnergy(int x, int y)
{
	double E_contrib = 0.0;
	//the modulo division ensures the periodic boundary conditions are met
	int right	= (x+1)%N;
	int left	= (x+N-1)%N;
	int down	= (y+1)%N;
	int up		= (y+N-1)%N;
	//add neighbour interaction energy to total energy
	E_contrib += spin[x][y]*spin[right][y];
	E_contrib += spin[x][y]*spin[x][down];
	E_contrib += spin[x][y]*spin[left][y];
	E_contrib += spin[x][y]*spin[x][up];

	return E_contrib *= -J;
}

double calcDeltaEnergy(int x, int y)
{
	return 2.0*calcSiteEnergy(x,y);
}

/*
*	MAGNETISATION FUNCTIONS
*/

double calcTotalMagnetisation()
{
	double M = 0.0;
	int right, down, left, up;
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			M += spin[x][y];
		}
	}
	M *= pow(N, -2.0);
	//multiply by factor at front of equation at end because I am an efficient coder LOLJK
	cout << "Magnetisation: " << M << "\n";
	return M;
}

/*
* HELPER FUNCTIONS
*/
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