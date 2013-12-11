// C/C++ header files
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <string>
// GNU Scientific Library header files
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
// my own header files
#include "config.h"
#include "ASCII.h"
#include "dimensions.cpp"

#define simulatedTime 100.0												//simulated time (seconds)
#define g 9.81																//acceleration due to gravity
#define pi 4.0*atan(1.0)													//mutha fuckin pi man
#define k_b 1.3806488e-23												//boltzmann's constant
#define mu 9.27400968e-24 												//bohr magneton

using namespace std;

/**** SHAMELESS GLOBAL VARIABLES ******/
int h = 0, c = 1;
int mesh[2][N][N];																//mesh[h] hot start, mesh[c] cold start

double total_E[2], total_M[2];
double E_av[2], E_s_av[2];												//average energy, average square energy
double S_av[2], S_s_av[2];												//average spin, average square spin

double previousEnergies[2][2][numberPrevEs];			//used to record moving averages of energy
double J = 1.0;												//J/(k_b*T) -> 1/(k_b*T) = J*beta
float mu_B = 0.5*J;

/**** SHAMELESS GLOBAL FILESTREAMS ****/
ofstream sys_props_file("data/sys_props.csv");
int seeds[10] = {116426264,1731462758, 1831960109, 851277867, 2075181264, 3079435518, 2241738047, 39641460, 4284274333, 1658052039};

/**** FUNCTION PROTOTYPES *************/

void simulateToEquilibrium(double);
void simulateXTimes(double, int);
void calculateSystemProperties();
void findEquilibrium(double, int);
bool atEquilibrium(int, double, int);

//spin functions
void initialiseSpins();
void initColdSpins();
void initHotSpins();
int randomSpin();
void flipSpin(int, int, int);

int calcTotalSpin(int);
int calcAverageSpin(int);
int calcTotalSpinSquared(int);
int calcAverageSpinSquared(int);

//random number functions
gsl_rng* setupUniformRNG(int i);
int randInt(int);

//energy functions
double calcAverageEnergy(int);
double calcTotalEnergy(int);
double calcPartialSiteEnergy(int, int, int);
double calcDeltaEnergy(int, int, int);

//magnetisation functions
double calcAverageMagnetisation(int);
double calcTotalMagnetisation(int);
double calcdM(int, int, int);

//
double calcSpecificHeatCapacity(double, int);
double calcMagneticSusceptibility(double, int);

//output to file/screen
void initOutputFile();
void outputSystemPropertiesToFile(double);

//helper functions
void updateProgress(double); 	//expects a fraction of completion
void done();


//for gsl rng help see here:
//http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-initialization.html#Random-number-generator-initialization
gsl_rng * r_uni;						//default rng instance


int main()
{
	r_uni = setupUniformRNG(2);
	initOutputFile();
	initialiseSpins();
	float beta;
	for (int i = 0; i <= int(beta_max/beta_step); ++i)
	{
		beta = i*beta_step;
		simulateToEquilibrium(beta);
		simulateXTimes(beta, N*N);
		calculateSystemProperties();
		outputSystemPropertiesToFile(beta);
		updateProgress(float(i*beta_step/beta_max));
	}
	done();
}

void simulateToEquilibrium(double beta)
{
	bool equilibrium[2] = {false, false};
	int i=0, i_outputted=0;									//the spin coordinates particular loop iteration

	total_E[h] = calcTotalEnergy(h);				//total energy of system
	total_M[h] = calcTotalMagnetisation(h);	//total magnetisation of system
	total_E[c] = calcTotalEnergy(c);				//total energy of system
	total_M[c] = calcTotalMagnetisation(c);	//total magnetisation of system

	//find equilibrium from hot start
	while(!equilibrium[h])
	{
		findEquilibrium(beta, h);
		equilibrium[h] = atEquilibrium(i, total_E[h], h);
		i++;
	}	//end of while loop for hot start

	//find equilibrium from hot start
	while(!equilibrium[c])
	{
		findEquilibrium(beta, c);
		equilibrium[c] = atEquilibrium(i, total_E[c], c);
		i++;
	}	//end of while loop for hot start

	//messages etc.
	// cout << "For beta = " << beta << ", equilibrium was found after " << i << " iterations.\n";
}

void simulateXTimes(double beta, int count)
{
	//once equilibrium has been reached we run the simulation
	//a few more times to [[[WHY?]]]
	for (int i = 0; i < count; ++i)
	{
		findEquilibrium(beta, h);
		findEquilibrium(beta, c);
	}
}

void calculateSystemProperties()
{
	total_E[h]	= calcTotalEnergy(h);
	total_M[h]	= calcTotalMagnetisation(h);
	E_av[h] 		= calcAverageEnergy(h);
	E_s_av[h] 	= 0.0;
	S_av[h] 		= calcAverageSpin(h);
	S_s_av[h]		= calcAverageSpinSquared(h);

	total_E[c]	= calcTotalEnergy(c);
	total_M[c]	= calcTotalMagnetisation(c);
	E_av[c] 		= calcAverageEnergy(c);
	E_s_av[c] 	= 0.0;
	S_av[c] 		= calcAverageSpin(c);
	S_s_av[c]		= calcAverageSpinSquared(c);
}

void findEquilibrium(double beta, int t)
{
	int x = randInt(N);
	int y = randInt(N);
	double dE = calcDeltaEnergy(x,y,t);
	if(dE < 0)
	{
		flipSpin(x,y,t);

		total_E[t] += dE;
		total_M[t] += calcdM(x,y,t);
	}
	else
	{
		double r = gsl_rng_uniform(r_uni); //random number 0 < r < 1
		if(r < exp(-dE*beta) )
		{
			flipSpin(x,y,t);
			
			total_E[t] += dE;
			total_M[t] += calcdM(x,y,t);
		}
	}
}

//returns true if mesh is at equilibrium, false otherwise
bool atEquilibrium(int i, double E, int t)
{
	//double pc_diff_max = 0.0001;

	if(i%(2*numberPrevEs) < numberPrevEs)
	{
		previousEnergies[t][0][i%numberPrevEs] = E;
	} else
	{
		previousEnergies[t][1][i%numberPrevEs] = E;
	}
	if( (i%numberPrevEs == 0) && (i>numberPrevEs) ) 
	{
		double	average1 = 0.0, average2 = 0.0;
		for (int i = 0; i < numberPrevEs; ++i)
		{
			average1 += previousEnergies[t][0][i];
			average2 += previousEnergies[t][1][i];
		}
		average1 /= numberPrevEs;
		average2 /= numberPrevEs;
		double pc_diff = abs((average1 - average2)/average1);							//percentage difference in average energies
		if(pc_diff < pc_diff_max)
		{
			return true;
		} else
		{
			return false;
		}
	} else
	{
		return false;
	}
}

gsl_rng* setupUniformRNG(int i)
{
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r,seeds[i]);
	return r;
}

int randInt(int max)
{
	return floor( gsl_rng_uniform(r_uni)*max );
}

/*
* SPIN MANIPULATION
*/

void initialiseSpins()
{
	initColdSpins();
	initHotSpins();
}

void initColdSpins()
{
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			mesh[c][x][y] = 1;
		}
	}
}

void initHotSpins()
{
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			mesh[h][x][y] = randomSpin();
		}
	}
}

int randomSpin()
{
	return gsl_rng_uniform(r_uni) > 0.5 ? 1 : -1 ;
}

void flipSpin(int x, int y, int t)
{
	mesh[t][x][y] = -1*mesh[t][x][y];
}

int calcTotalSpin(int t)
{
	int s = 0;
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			s += mesh[t][x][y];
		}
	}
	return s;
}

int calcAverageSpin(int t)
{
	return calcTotalSpin(t)/pow(N, 2.0);
}

int calcTotalSpinSquared(int t)
{
	int s = 0;
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			s += pow(mesh[t][x][y], 2.0);
		}
	}
	return s;
}

int calcAverageSpinSquared(int t)
{
	return calcTotalSpinSquared(t)/pow(N, 2.0);
}

/*
*	ENERGY FUNCTIONS
*/
double calcAverageEnergy(int t)
{
	return calcTotalEnergy(t)/pow(N, 2.0);
}

double calcTotalEnergy(int t)
{
	double E_1 = 0.0,  E_2 = 0.0;
	//int right, down, left, up;
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			E_1 += calcPartialSiteEnergy(x,y,t);
			E_2 += mesh[t][x][y];
		}
	}
	//multiply by factor at front of equation at end because I am an efficient coder LOLJK
	return -0.5*J*E_1 - mu_B*E_2;
}

double calcPartialSiteEnergy(int x, int y, int t)
{
	double E_contrib = 0.0;
	//the modulo division ensures the periodic boundary conditions are met
	int right	= (x+1)%N;
	int left	= (x+N-1)%N;
	int down	= (y+1)%N;
	int up		= (y+N-1)%N;
	//add neighbour interaction energy to total energy
	E_contrib += mesh[t][x][y]*mesh[t][right][y];
	E_contrib += mesh[t][x][y]*mesh[t][x][down];
	E_contrib += mesh[t][x][y]*mesh[t][left][y];
	E_contrib += mesh[t][x][y]*mesh[t][x][up];

	return E_contrib *= J;
}

double calcDeltaEnergy(int x, int y, int t)
{
	return 2.0*calcPartialSiteEnergy(x,y,t);
}

/*
*	MAGNETISATION FUNCTIONS
*/

double calcAverageMagnetisation(int t)
{
	return calcTotalMagnetisation(t)/pow(N, 2.0);
}

double calcTotalMagnetisation(int t)
{
	double M = 0.0;
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			M += mesh[t][x][y];
		}
	}
	M *= pow(N, -2.0);
	//multiply by factor at front of equation at end because I am an efficient coder LOLJK
	//cout << "Magnetisation: " << M << "\n";
	return M;
}

double calcdM(int x, int y, int t)
{
	return 2*mesh[t][x][y]*pow(N, -2.0);
}

/*
*
*/
double calcSpecificHeatCapacity(double beta, int t)
{
	//this one I think is right
	//return pow(N, -2.0) * ( (k_b*pow(beta,2.0)) / pow(J,2.0) ) * ( E_s_av[t] - pow(E_av[t], 2.0)) ;
	return pow(N, -2.0) * ( pow(beta,2.0) / (k_b*pow(J,2.0)) ) * ( E_s_av[t] - pow(E_av[t], 2.0)) ;
}

double calcMagneticSusceptibility(double beta, int t)
{
	return pow(N, -2.0)*(beta) * ( S_s_av[t] - pow(S_av[t], 2.0));
}


/*
*	OUTPUT TO FILE/SCREEN
*/
void initOutputFile()
{
	sys_props_file << "beta";
	if(outputE) sys_props_file << ",E_av (cold),E_av (hot)";
	if(outputM) sys_props_file << ",M_av (cold),M_av (hot)";
	if(outputSHC) sys_props_file << ",SHC (cold),SHC (hot)";
	if(outputMS) sys_props_file << ",MS (cold),MS (hot)";
	sys_props_file << "\n";
}

void outputSystemPropertiesToFile(double beta)
{
	sys_props_file	<< beta;
	if(outputE) 	sys_props_file << "," << total_E[c]/pow(N, 2.0) << "," << total_E[h]/pow(N, 2.0);
	if(outputM) 	sys_props_file << "," << total_M[c]/pow(N, 2.0) << "," << total_M[h]/pow(N, 2.0);
	if(outputSHC) sys_props_file << "," << calcSpecificHeatCapacity(beta, c) << "," << calcSpecificHeatCapacity(beta, h);
	if(outputMS) 	sys_props_file << "," << calcMagneticSusceptibility(beta, c) << "," << calcMagneticSusceptibility(beta, h);
	sys_props_file << "\n";
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

	for(int c = 0; c < progressBarLength; ++c)
	{
		progressBar << "#";
	}
	for(int c = 0; c < width-progressBarLength; ++c)
	{
		progressBar << " ";
	}
	//std::cout << "\rRunning: " << progress*100 << "\%\n" << std::flush;
	std::cout << "\r" << progressBar.str() << std::flush;
}

//displays done message to terminal
void done()
{
	std::cout.precision(3);
	std::cout << std::fixed << "\r\n\n"
						<< ASCII_done
						<<	"Program Report:\n"
						<< "---------------\n"
						<< "Mesh dimensions: " << N << "x" << N << "\n"
						<< "Beta range explored: 0 - " << beta_max << " in steps of " << beta_step << "\n"
						<<	"J: " << J << "\tmu_B: " << mu_B << "\n"
						<< "==========================================\n"
						<< "Property\t\t\tOutputted\n"
						<< "------------------------------------------\n"
						<< "Energy\t\t\t\t" << (outputE ? "YES\n" : "NO\n")
						<< "Magnetisation\t\t\t" << (outputM ? "YES\n" : "NO\n")
						<< "Specific Heat Capacity\t\t" << (outputSHC ? "YES\n" : "NO\n")
						<< "Magnetic Susceptibility\t\t" << (outputMS ? "YES\n" : "NO\n")
						<< "==========================================\n";
}