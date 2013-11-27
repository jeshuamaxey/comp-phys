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
#define numberPrevEs 50														//number of previous energies we keep track of to test for equilibrium
//#define beta 0																		// J/(k_b*T) -> 1/(k_b*T) = J*beta
#define beta_step 0.001
#define beta_max 1.0

/**** PROGRAM CONTROL SETTINGS ********/
#define coldStart false														//use for starting at low T

using namespace std;

/**** SHAMELESS GLOBAL VARIABLES ******/
double total_E, total_M;

/**** SHAMELESS GLOBAL FILESTREAMS ****/
ofstream sys_props_file_c("data/sys_props_cold.csv");
ofstream sys_props_file_h("data/sys_props_hot.csv");
ofstream json("data/mesh.json");

/**** FUNCTION PROTOTYPES *************/

void simulateToEquilibrium(double, double);
void simulateSomeMore(double);
void calculateSystemProperties();
void findEquilibrium(double);
bool atEquilibrium(int, double);

//spin functions
void initialiseSpins();
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
double calcdM(int, int);

//output to file/screen
void outputSystemPropertiesToFile(double);
void outputSpinsToFile();
void initJsonFile(double beta);
void updateJsonFile(int);
void endJsonFile(int);

//helper functions
void updateProgress(double);
void done();

//spin mesh
int spin[N][N];
double previousEnergies[2][numberPrevEs];
double J = 1.0;

//for gsl rng help see here:
//http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-initialization.html#Random-number-generator-initialization
gsl_rng * r_uni;						//default rng instance

int main()
{
	cout << "\nPROJECT-B.CPP\n=============\n\n";
	sys_props_file_c << "beta,M (cold),E (cold)\n";
	sys_props_file_h << "beta,M (hot),E (hot\n";
	double initial_temp = 0.0, beta;

	r_uni = setupUniformRNG();
	initialiseSpins();
	for (int i = 0; i <= int(beta_max/beta_step); ++i)
	{
		beta = i*beta_step;
		simulateToEquilibrium(initial_temp, beta);
		simulateSomeMore(beta);
		calculateSystemProperties();
		outputSystemPropertiesToFile(beta);
		//updateProgress(float(i*beta_max/beta_step));
	}
	outputSpinsToFile();
	done();
}

void simulateToEquilibrium(double initial_temp, double beta)
{
	//
	initJsonFile(beta);
	bool equilibrium = false;
	int i=0, i_outputted=0;							//the spin coordinates particular loop iteration

	total_E = calcTotalEnergy();				//total energy of system
	total_M = calcTotalMagnetisation();	//total magnetisation of system

	//find equilibrium
	while(!equilibrium)
	{
		findEquilibrium(beta);
		equilibrium = atEquilibrium(i, total_E);
		/*
		if(i%10 == 0)
		{
			updateJsonFile(i);
			i_outputted++;
		}
		*/
		i++;
	}	//end of while loop

	
	//messages etc.
	// cout << "equilibrium found after " << i << " iterations.\n";
	// cout << i_outputted << " iterations outputted to JSON file.\n";
	endJsonFile(i_outputted);
}

void simulateSomeMore(double beta)
{
	//once equilibrium has been reached we run the simulation
	//a few more times to [[[WHY?]]]
	for (int i = 0; i < N*N; ++i)
	{
		findEquilibrium(beta);
	}
}

void calculateSystemProperties()
{
	total_E = calcTotalEnergy();
	total_M = calcTotalMagnetisation();
	// E_av, E_s_av,						//average energy, average square energy
	// S_av, S_s,av;						//average spin, average square spin
}

void outputSystemPropertiesToFile(double beta)
{
	if(coldStart) sys_props_file_c << beta << "," << total_M << "," << total_E << "\n";
	else 					sys_props_file_h << beta << "," << total_M << "," << total_E << "\n";
}

void findEquilibrium(double beta)
{
	int x = randInt(N);
	int y = randInt(N);
	double dE = calcDeltaEnergy(x,y);
	if(dE < 0)
	{
		flipSpin(x,y);

		total_E += dE;
		total_M += calcdM(x,y);
	}
	else
	{
		double r = gsl_rng_uniform(r_uni); //random number 0 < r < 1
		if(r < exp(-dE*beta) )
		{
			flipSpin(x,y);
			
			total_E += dE;
			total_M += calcdM(x,y);
		}
	}
}

//returns true if mesh is at equilibrium, false otherwise
bool atEquilibrium(int i, double E)
{
	double pc_diff_max = 0.001;
	if(i%(2*numberPrevEs) < numberPrevEs)
	{
		previousEnergies[0][i%numberPrevEs] = E;
	} else
	{
		previousEnergies[1][i%numberPrevEs] = E;
	}
	if( (i%numberPrevEs == 0) && (i>numberPrevEs) ) 
	{
		double	average1 = 0.0, average2 = 0.0;
		for (int c = 0; c < numberPrevEs; ++c)
		{
			average1 += previousEnergies[0][c];
			average2 += previousEnergies[1][c];
		}
		average1 /= numberPrevEs;
		average2 /= numberPrevEs;
		double pc_diff = abs((average1 - average2)/average1);							//percentage difference in average energies
		if(pc_diff < pc_diff_max)
		{
			/*
			cout 	<< "Equilibrium conditions:\n-----------------------\n";
			cout 	<< "av_1\t\tav_2\t\tpc_diff\t\tpc_diff_max\tno. prev E.s\n";
			cout 	<< "-----\t\t-----\t\t-----\t\t-----\t\t-----\t\t\n";
			cout 	<< average1 << "\t\t"
						<< average2 << "\t\t"
						<< pc_diff*100 	<< "%\t\t"
						<< pc_diff_max*100 << "%\t\t"
						<< numberPrevEs << "\n\n";
			*/
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
			E += 0.5*J*calcSiteEnergy(x,y);
		}
	}
	//multiply by factor at front of equation at end because I am an efficient coder LOLJK
	//cout << "Energy: " << E << "\n";
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
	//cout << "Magnetisation: " << M << "\n";
	return M;
}

double calcdM(int x, int y)
{
	return 2*spin[x][y]*pow(N, -2.0);
}

/*
*	OUTPUT TO FILE/SCREEN
*/
void initJsonFile(double beta)
{
	json << "{\n\"mesh\" : [\n";
}
void updateJsonFile(int i)
{
	stringstream mesh;
	//mesh << "[";
	for (int x = 0; x < N; ++x)
	{
		mesh << "[";
		for (int y = 0; y < N; ++y)
		{
			mesh << spin[x][y];
			if(y!=N-1) mesh << ",";
		}
		if(x!=N-1) mesh << "],\n\t";
		else mesh << "]\n";
	}
	if(i == 0) json << "[\n\t"<< mesh.str() << "]\n";
	else json << ",[\n\t"<< mesh.str() << "]\n";
}

void endJsonFile(int iterations) 
{
	json << "],\n\"config\" : {\n";
	json << "\"iterations\" : " << iterations << ",\n";
	json << "\"meshSize\" : " << N << "\n}\n}";
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
	std::cout << "\r" << progressBar.str() << std::flush;
}

//displays done message to terminal
void done()
{
	std::cout << std::fixed << "\r"<< "\n100\% - all done. Simulated time elapsed: " << simulatedTime << "s\n" << std::flush;
}