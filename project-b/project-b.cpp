// C/C++ header files
#include <cstdlib>
#include <stdint.h>
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
#include "tracktime.cpp"

#define simulatedTime 100.0												//simulated time (seconds)
#define g 9.81																		//acceleration due to gravity
#define pi 4.0*atan(1.0)													//mutha fuckin pi man
#define k_b 1.3806488e-23													//boltzmann's constant
#define mu 9.27400968e-24 												//bohr magneton

using namespace std;

/**** SHAMELESS GLOBAL VARIABLES ******/
int h = 0, c = 1;
int mesh[2][N][N];																//mesh[h] hot start, mesh[c] cold start

double 	totalMeshEnergiesPastEquilibrium[2][simulationsPastEquilibrium],
				totalMeshSpinPastEquilibrium[2][simulationsPastEquilibrium],
				totalMeshMagnetisationPastEquilibrium[2][simulationsPastEquilibrium];

double total_mesh_E[2], total_mesh_M[2];											//total energy, total magnetisation (of micro state)
double E_av[numberOfSeeds][2], E_s_av[numberOfSeeds][2];			//average energy, average square energy
double S_av[numberOfSeeds][2], S_s_av[numberOfSeeds][2];			//average spin, average square spin
double magneticSusceptibility[numberOfSeeds][2];							//
double specificHeatCapacity[numberOfSeeds][2];								//

double previousEnergies[2][2][numberPrevEs];			//used to record moving averages of energy

//char processName[3][64] = {'Simulating zero B','Simulating non-zero B','Simulating varying B'};

int seeds[numberOfSeeds] = {
		116426264,
		731462758,
		1831960109,
		851277867,
		2075181264,
		3079435518,
		2241738047,
		39641460,
		4284274333,
		1658052039 };

double startTime, endTime;												//keep track of how long code took to run

/**** SHAMELESS GLOBAL FILESTREAMS ****/
ofstream sys_props_file("data/sys_props.csv");
ofstream debug("data/debug.csv");

/**** FUNCTION PROTOTYPES *************/
void initOutputFile(ostream&);

void runSimulation(double, double, ostream&, int, int);

void simulateZeroB(ostream&, int);
void simulateNonZeroB(ostream&, int);
void simulateVaryingB(ostream&, int);

void simulateToEquilibrium(double, double);
void simulatePastEquilibrium(double, double);
void calculateSystemProperties(int, double);
void metropolisSpinFlip(double, int);
bool atEquilibrium(int, double, int);

//spin functions
void initialiseSpins();
void initColdSpins();
void initHotSpins();
int randomSpin();
void flipSpin(int, int, int);

int calcTotalSpin(int);
double calcAverageSpin(int);
double calcAverageSpinSquared(int);

//random number functions
gsl_rng* setupUniformRNG(int i);
int randInt(int);

//energy functions
double calcAverageEnergy(int);
double calcAverageEnergySquared(int);
double calcTotalMicroEnergy(int, double);
double calcPartialSiteEnergy(int, int, int);
double calcDeltaEnergy(int, int, int);
double calcMeanAverageEnergy(int);

//magnetisation functions
double calcTotalMicroMagnetisation(int);
double calcdM(int, int, int);
double calcTotalMacroMagnetisation(double, int);

//specific heat capacity functions
double calcSpecificHeatCapacity(double, int, int);
double calcMeanSpecificHeatCapacity(double, int);

//magnetic susceptibility functions
double calcMagneticSusceptibility(double, int, int);
double calcMeanMagneticSusceptibility(double, int);

//output to file/screen
void outputSystemPropertiesToFile(double, double, ostream&, int);

//helper functions
void updateProgress(double, int); 	//expects a fraction of completion
void done();


//for gsl rng help see here:
//http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-initialization.html#Random-number-generator-initialization
gsl_rng * r_uni;						//default rng instance


int main()
{
	//log start time
	startTime = GetTimeMs64();
	//create filestreams
	ofstream zeroBFile("data/zeroB.csv");
	ofstream nonZeroBFile("data/nonZeroB.csv");
	ofstream varyingZeroBFile("data/varyingZeroB.csv");

	//absolutely necessary debugging tool
	cout << ASCII_batman;

	initOutputFile(zeroBFile);
	initOutputFile(nonZeroBFile);
	initOutputFile(varyingZeroBFile);
	//loop that iterates through rng seeds
	for (int a = 0; a < numberOfSeeds; ++a)
	{
		//reset random number generator
		r_uni = setupUniformRNG(a);
		//align spins
		initialiseSpins();
		//run appropriate simulations
		if(simulateForZeroB) simulateZeroB(zeroBFile, a);
		if(simulateForNonZeroB) simulateNonZeroB(nonZeroBFile, a);
		if(simulateForVaryingB) simulateVaryingB(varyingZeroBFile, a);
	}

	//job up!
	done();
}

void simulateZeroB(ostream& outputFile, int a)
{
	//Zero B Field - varying beta
	//initOutputFile(outputFile);
	//fix mu_B
	float beta =0, mu_B = 0;
	for (int i = 0; i <= int(beta_max/beta_step); ++i)
	{
		//set beta
		beta = beta_min+(i*beta_step);
		//run simulation under these conditions
		runSimulation(beta, mu_B, outputFile, a, 0);
		//
		updateProgress(float(i)/(beta_max/beta_step), a);
	}
}

void simulateNonZeroB(ostream& outputFile, int a)
{
	//Constant (non-zero) B Field - varying beta
	//initOutputFile(outputFile);
	//fix mu_B
	float beta =0, mu_B = 0.5*J;
	for (int i = 0; i <= int(beta_max/beta_step); ++i)
	{
		//set beta
		beta = beta_min+(i*beta_step);
		//run simulation under these conditions
		runSimulation(beta, mu_B, outputFile, a, 1);
		//
		updateProgress(float(i)/(beta_max/beta_step), a);
	}
}

void simulateVaryingB(ostream& outputFile, int a)
{
	//Constant beta (~T_crit ) - varying B field
	//initOutputFile(outputFile);
	//fix beta
	float beta = 0.25, mu_B = 0;
	for (int i = 0; i <= int(mu_B_max/mu_B_step); ++i) {
		//set magnetic field strength
		mu_B = mu_B_min+(i*mu_B_step);
		//run simulation under these conditions
		runSimulation(beta, mu_B, outputFile, a, 2);
		//
		updateProgress(float(i)/(mu_B_max/mu_B_step), a);
	}
}

void runSimulation(double beta, double mu_B, ostream& outputFile, int a, int nameIndex)
{
	simulateToEquilibrium(beta, mu_B);
	simulatePastEquilibrium(beta, mu_B);
	calculateSystemProperties(a, beta);
	//on last iteration of loop output to file
	if(a==numberOfSeeds-1) outputSystemPropertiesToFile(beta, mu_B, outputFile, a);
}

void simulateToEquilibrium(double beta, double mu_B)
{
	bool equilibrium[2] = {false, false};
	int i=0, i_outputted=0;																	//the spin coordinates particular loop iteration

	total_mesh_E[h] = calcTotalMicroEnergy(h, mu_B);				//total energy of system
	total_mesh_M[h] = calcTotalMicroMagnetisation(h);				//total magnetisation of system
	total_mesh_E[c] = calcTotalMicroEnergy(c, mu_B);				//total energy of system
	total_mesh_M[c] = calcTotalMicroMagnetisation(c);				//total magnetisation of system

	//find equilibrium from hot start
	while(!equilibrium[h])
	{
		metropolisSpinFlip(beta, h);
		equilibrium[h] = atEquilibrium(i, total_mesh_E[h], h);
		i++;
	}	//end of while loop for hot start

	i=0;
	//find equilibrium from cold start
	while(!equilibrium[c])
	{
		metropolisSpinFlip(beta, c);
		equilibrium[c] = atEquilibrium(i, total_mesh_E[c], c);
		i++;
	}	//end of while loop for cold start
}

void simulatePastEquilibrium(double beta, double mu_B)
{
	//once equilibrium has been reached we run the simulation
	//a few more times to [[[WHY?]]]
	//arrays for
	for(int i = 0; i < simulationsPastEquilibrium; ++i)
	{
		//
		for(int j = 0; j < N*N; ++j)
		{
			metropolisSpinFlip(beta, h);
			metropolisSpinFlip(beta, c);
		}
		totalMeshEnergiesPastEquilibrium[h][i] = calcTotalMicroEnergy(h, mu_B);
		totalMeshSpinPastEquilibrium[h][i] = calcTotalSpin(h);
		totalMeshMagnetisationPastEquilibrium[h][i] = calcTotalMicroMagnetisation(h);

		totalMeshEnergiesPastEquilibrium[c][i] = calcTotalMicroEnergy(c, mu_B);
		totalMeshSpinPastEquilibrium[c][i] = calcTotalSpin(c);
		totalMeshMagnetisationPastEquilibrium[c][i] = calcTotalMicroMagnetisation(c);
	}
}

void calculateSystemProperties(int a, double beta)
{
	E_av[a][h] 		= calcAverageEnergy(h);
	E_s_av[a][h] 	= calcAverageEnergySquared(h);
	S_av[a][h] 		= calcAverageSpin(h);
	S_s_av[a][h]	= calcAverageSpinSquared(h);

	E_av[a][c] 		= calcAverageEnergy(c);
	E_s_av[a][c] 	= calcAverageEnergySquared(c);
	S_av[a][c] 		= calcAverageSpin(c);
	S_s_av[a][c]	= calcAverageSpinSquared(c);

	//mag sus
	magneticSusceptibility[a][h] = calcMagneticSusceptibility(beta, h, a);
	magneticSusceptibility[a][c] = calcMagneticSusceptibility(beta, c, a);
	//shc
	specificHeatCapacity[a][h] = calcSpecificHeatCapacity(beta, h, a);
	specificHeatCapacity[a][c] = calcSpecificHeatCapacity(beta, c, a);
}

void metropolisSpinFlip(double beta, int t)
{
	int x = randInt(N);
	int y = randInt(N);
	double dE = calcDeltaEnergy(x,y,t);
	if(dE < 0)
	{
		flipSpin(x,y,t);

		total_mesh_E[t] += dE;
		total_mesh_M[t] += calcdM(x,y,t);
	}
	else
	{
		double r = gsl_rng_uniform(r_uni); //random number 0 < r < 1
		if(r < exp(-dE*beta) )
		{
			flipSpin(x,y,t);
			
			total_mesh_E[t] += dE;
			total_mesh_M[t] += calcdM(x,y,t);
		}
	}
}

//returns true if mesh is at equilibrium, false otherwise
bool atEquilibrium(int i, double E, int t)
{
	int x = i%(4*numberPrevEs);

	// 50% chance we need to ignore so do shortcircuiting costly evaluation first
	if( ( (0 <= x) && (x< 3*N) ) || (6*N <= x && x< 9*N) )	//	0 < x < 3N OR 6N < x < 9N
	{
		//ignore this value
		return false;
	}
	//is 3N < x < 6N ?
	//store in first moving average array
	else if (3*N <= x && x< 6*N)
	{
		previousEnergies[t][0][i%numberPrevEs] = E;
	}
	//is 9N < x < 12N
	//store in second moving average array
	else if (9*N <= x && x< 12*N)
	{
		previousEnergies[t][1][i%numberPrevEs] = E;
	}
	if( (i%(4*numberPrevEs -1) == 0) && (i>0) ) 
	{
		//calculate new moving averages
		double	average1 = 0.0, average2 = 0.0;
		for (int i = 0; i < numberPrevEs; ++i)
		{
			average1 += previousEnergies[t][0][i];
			average2 += previousEnergies[t][1][i];
		}
		average1 /= numberPrevEs;
		average2 /= numberPrevEs;
		//compare against pc_diff_max which characterises the equilibrium condition
		double pc_diff = abs((average1 - average2)/average1);			//percentage difference in average energies
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
	//program crashes here
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

double calcAverageSpin(int t)
{
	double s = 0.0;
	for (int i = 0; i < simulationsPastEquilibrium; ++i)
	{
		s += totalMeshSpinPastEquilibrium[t][i];
	}
	return s/float(simulationsPastEquilibrium);
}

double calcAverageSpinSquared(int t)
{
	int s = 0;
	for (int i = 0; i < simulationsPastEquilibrium; ++i)
	{
		s += pow(totalMeshSpinPastEquilibrium[t][i], 2.0);
	}
	return s/float(simulationsPastEquilibrium);
}

/*
*	ENERGY FUNCTIONS
*/
double calcAverageEnergy(int t)
{
	double E = 0.0;
	for (int i = 0; i < simulationsPastEquilibrium; ++i)
	{
		E += totalMeshEnergiesPastEquilibrium[t][i];
	}
	return E/float(simulationsPastEquilibrium);
}

double calcAverageEnergySquared(int t)
{
	double E = 0.0;
	for (int i = 0; i < simulationsPastEquilibrium; ++i)
	{
		E += pow(totalMeshEnergiesPastEquilibrium[t][i], 2.0);
	}
	return E/float(simulationsPastEquilibrium);
}

double calcTotalMicroEnergy(int t, double mu_B)
{
	double E_1 = 0.0,  E_2 = 0.0;
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

double calcMeanAverageEnergy(int t)
{
	double E = 0.0;
	for (int i = 0; i < numberOfSeeds ; ++i)
	{
		E += E_av[i][t];
	}
	return E/float(numberOfSeeds);
}

/*
*	MAGNETISATION FUNCTIONS
*/

double calcTotalMicroMagnetisation(int t)
{
	double M = 0.0;
	for (int x = 0; x < N; ++x)
	{
		for (int y = 0; y < N; ++y)
		{
			M += mesh[t][x][y];
		}
	}
	//multiply by factor at front of equation at end because I am an efficient coder LOLJK
	M *= pow(N, -2.0);
	return M;
}

double calcdM(int x, int y, int t)
{
	return 2*mesh[t][x][y]*pow(N, -2.0);
}

double calcTotalMacroMagnetisation(double beta, int t)
{
	double M = 0.0;
	for (int i = 0; i < simulationsPastEquilibrium; ++i)
	{
		M += abs(totalMeshMagnetisationPastEquilibrium[t][i]);
	}
	return M/(N*N);
}

/*
* Specific heat capacity formulae
*/
double calcSpecificHeatCapacity(double beta, int t, int a)
{
	return pow(N, -2.0) * ( (k_b*pow(beta,2.0)) / pow(J,2.0) ) * ( E_s_av[a][t] - pow(E_av[a][t], 2.0));
}

double calcMeanSpecificHeatCapacity(double beta, int t)
{
	double shc = 0.0;
	for (int i = 0; i < numberOfSeeds; ++i)
	{
		shc += specificHeatCapacity[i][t];
	}
	return shc/float(numberOfSeeds);
}

/*
* magnetic susceptibility formulae
*/

double calcMagneticSusceptibility(double beta, int t, int a)
{
	return pow(N, -2.0)*(beta) * ( S_s_av[a][t] - pow(S_av[a][t], 2.0));
}

double calcMeanMagneticSusceptibility(double beta, int t)
{
	float ms = 0.0;
	for (int i = 0; i < numberOfSeeds; ++i)
	{
		ms += magneticSusceptibility[i][t];
	}
	return ms/float(numberOfSeeds);
}


/*
*	OUTPUT TO FILE/SCREEN
*/

void initOutputFile(ostream& outputFile)
{
	if(outputOneOverBeta) outputFile << "T";
	else 									outputFile << "beta";
	if(outputE) outputFile << ",E_av (cold),E_av (hot)";
	if(outputM) outputFile << ",M (cold),M (hot)";
	if(outputSHC) outputFile << ",SHC (cold),SHC (hot)";
	if(outputMS) outputFile << ",MS (cold),MS (hot)";
	outputFile << "\n";
}

void outputSystemPropertiesToFile(double beta, double mu_B, ostream& outputFile, int a)
{
	if(outputOneOverBeta) outputFile	<< 1/beta;
	else 									outputFile	<< beta;
	if(outputE) 	outputFile << "," << calcMeanAverageEnergy(c) << "," << calcMeanAverageEnergy(h);
	if(outputM) 	outputFile << "," << calcTotalMacroMagnetisation(beta, c) << "," << calcTotalMacroMagnetisation(beta, h);
	if(outputSHC) outputFile << "," << calcMeanSpecificHeatCapacity(beta, c) << "," << calcMeanSpecificHeatCapacity(beta, h);
	if(outputMS) 	outputFile << "," << calcMeanMagneticSusceptibility(beta, c) << "," << calcMeanMagneticSusceptibility(beta, h);
	outputFile << "\n";
}

/*
* HELPER FUNCTIONS
*/
//updates the progress display on commandline
void updateProgress(double progress, int a)
{
	int width = getTerminalWidth();
	progress /= float(numberOfSeeds);
	progress += float(a)/float(numberOfSeeds);
	int progressBarLength = int(progress*width);
	stringstream progressBar;

	for(int c = 0; c < progressBarLength; ++c)
	{
		progressBar << a;
	}
	for(int c = 0; c < width-progressBarLength; ++c)
	{
		progressBar << " ";
	}
	std::cout << "\r" << progressBar.str() << std::flush;
}

//displays done message to terminal
void done()
{
	std::cout.precision(3);
	std::cout << std::fixed << "\r\n\n"
						<< ASCII_program_report
						<< "-------------------------------------------------------------\n"
						<< "Time to execute: " << (GetTimeMs64() - startTime)/1000 << " seconds\n"
						<< "Mesh dimensions: " << N << "x" << N << "\n"
						<< "Beta range explored: "<< beta_min<< " - " << beta_max << " in steps of " << beta_step << "\n"
						<< "mu_B range explored: "<< mu_B_min<< " - " << mu_B_max << " in steps of " << mu_B_step << "\n"
						<< "Beta outputted as temp: " << (outputOneOverBeta ? "OBVIOUSLY\n" : "YOU HAVING A LAUGH\n")
						<<	"J: " << J << "\n"
						<< "-------------------------------------------------------------\n"
						<< "=============================================================\n"
						<< "Simulation\t\t\t\t\tRun\n"
						<< "-------------------------------------------------------------\n"
						<< "Varying beta, no magnetic field\t\t\t" << (simulateForZeroB ? "YES\n" : "NO\n")
						<< "Varying beta, constant magnetic field\t\t" << (simulateForNonZeroB ? "YES\n" : "NO\n")
						<< "Constant beta, varying magnetic field\t\t" << (simulateForVaryingB ? "YES\n" : "NO\n")
						<< "=============================================================\n\n"
						<< "=============================================================\n"
						<< "Property\t\t\t\t\tOutputted\n"
						<< "-------------------------------------------------------------\n"
						<< "Energy\t\t\t\t\t\t" << (outputE ? "YES\n" : "NO\n")
						<< "Magnetisation\t\t\t\t\t" << (outputM ? "YES\n" : "NO\n")
						<< "Specific Heat Capacity\t\t\t\t" << (outputSHC ? "YES\n" : "NO\n")
						<< "Magnetic Susceptibility\t\t\t\t" << (outputMS ? "YES\n" : "NO\n")
						<< "=============================================================\n";
}