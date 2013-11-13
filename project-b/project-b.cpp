#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <string>
#include "dimensions.cpp"

#define simulatedTime 100.0												//simulated time (seconds)
#define g 9.81																		//acceleration due to gravity
#define pi atan(1.0)															//mutha fuckin pi man

#define N 10																			//there are NxN spins simulated

/****** PROGRAM CONTROL SETTINGS ******/
#define initSpinsAligned true											//use for starting at low T
#define initSpinsRandom false											//use for starting at high T

using namespace std;

/****** FUNCTION PROTOTYPES ***********/
void initialiseSpins();
void findEquilibrium(double);
void alignSpins();
void randomlyDistSpins();
int randIndex(int);

void updateProgress(double, char[]);
void done();

int spin[N][N];

int main()
{
	double initial_temp = 0.0;

	initialiseSpins();
	findEquilibrium(initial_temp);

	long int longtime = 100000;
	char str[64] = "pooing";

	for (int i = 0; i < longtime; ++i)
	{
		updateProgress((1.0*i/longtime), str);
	}
}

/* START ALL FUNCTIONS CALLED FROM main() */

void initialiseSpins()
{
	if(initSpinsAligned)
	{
		alignSpins();
	}
	else if(initSpinsRandom)
	{
		randomlyDistSpins();
	}
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

/* END ALL FUNCTIONS CALLED FROM main() */

void alignSpins()
{}

void randomlyDistSpins()
{}

int randIndex(int max)
{
	return 1;
}

//updates the progress display on commandline
void updateProgress(double progress, char *name)
{
	int width = getTerminalWidth();
	int progressBarLength = int(progress*width);
	//std::cout << progressBarLength;
	stringstream progressBar;

	for(int c=0;c<progressBarLength;c++)
	{
		progressBar << "#";
	}
	std::cout << "\r" << progressBar.str() << std::flush;
	//std::cout << "\r"<< "Currently running " << name << " - " << progress << "\t" << int(progress*100) << "%" << std::flush;
}

//displays done message to terminal
void done()
{
	std::cout << std::fixed << "\r"<< "100\% - all done. Simulated time elapsed: " << simulatedTime << "s\n" << std::flush;
}