#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <string>

#define simulatedTime 10.0											//simulated time (seconds)
#define h 0.01																	//step size
#define numberOfSteps  int(simulatedTime / h)		//used for sizing arrays
#define g 9.81																	//acceleration due to gravity
#define pi atan(1.0)														//mutha fuckin pi man

#define outputPositions false										//determine what's outputted by the program to file
#define outputEnergies true											//determine what's outputted by the program to file

using namespace std;

void initFile(ostream&, string);
string makePositionHeading(string);
string makeEnergyHeading(string);
void outputPositionToFile(ostream&, double, double);
void outputEnergyToFile(ostream&, double, double, double);

void setInitialValues();

void single_pendulum();

void updateEnergies(int);
void updateEuler(int);
void updateLeapfrog(int);
void updateRK4(int);

double calculateKineticEnergy(double m, double l, double theta, double w);
double calculatePotentialEnergy(double m, double l, double theta, double w);

void eulerStabilityTest();
void stabilityTest();
string isStable(double, double);
void calculateAnalyticalSolution();

//helper functions
void updateProgress(int);
void updateProgress(int, char[]);
void done();

int main()
{
	setInitialValues();
	single_pendulum();
	stabilityTest();
	eulerStabilityTest();
	done();
	return 0;
}

/********** GLOBAL VARIABLES COS FUCK EFFICIENCY ***********/

double t = 0;																//time
double l = 9.81;														//length of pendulem in metres
double m = 1.0;															//mass of pendulum in kg
double damping_constant = 0.0;													//damping coefficient
double beta = damping_constant/(m* sqrt( g*l ));				//matrix constant

//initial conditions
double initial_theta = 0.15;									//angle from vert (starting angle)
double initial_w = 0;												//rate of change of angle from vert (starting at rest)

//euler variable arrays
double euler_theta [numberOfSteps];					//stores all theta values
double euler_w [numberOfSteps];							//stores all w values
double euler_E [numberOfSteps];							//stores all E values
double euler_T [numberOfSteps];							//stores all T values
double euler_U [numberOfSteps];							//stores all U values

//leapfrog variable arrays
double leapfrog_theta [numberOfSteps];			//stores all theta values
double leapfrog_w [numberOfSteps];					//stores all w values
double leapfrog_E [numberOfSteps];					//stores all E values
double leapfrog_T [numberOfSteps];					//stores all T values
double leapfrog_U [numberOfSteps];					//stores all U values

//rk4 variable arrays
double rk4_theta [numberOfSteps];						//stores all theta values
double rk4_w [numberOfSteps];								//stores all w values
double rk4_E [numberOfSteps];								//stores all E values
double rk4_T [numberOfSteps];								//stores all T values
double rk4_U [numberOfSteps];								//stores all U values
double k_1, k_2, k_3, k_4;

//analytical solution & stabiliy testvars
double anal_theta[numberOfSteps];
double anal_w[numberOfSteps];
double err_theta[numberOfSteps];

/***************************************************************/

void initFile(ostream& file, string processName)
{
	//first column is always time
	file << "time" ;
	//position heading
	if(outputPositions)
	{
		file << makePositionHeading(processName);
	}
	//energy heading
	if(outputEnergies)
	{
		file << makeEnergyHeading(processName);
	}
	//end column headers
	file << "\n";
}

string makePositionHeading(string processName)
{
	stringstream string;
	string << "," << processName << "_theta," << processName << "_w" ;
	return string.str();
}

string makeEnergyHeading(string processName)
{
	stringstream string;
	string << "," << processName << "_T," << processName << "_U," << processName << "_E" ;
	return string.str();
}

void setInitialValues()
{
	//set initial values
	euler_theta[0] = initial_theta;
	euler_w[0] = initial_w;

	leapfrog_theta[0] = initial_theta;
	leapfrog_w[0] = initial_w;

	rk4_theta[0] = initial_theta;
	rk4_w[0] = initial_w;
}

/*
                                      
          88                         88                                                   
          ""                         88                                                   
                                     88                                                   
,adPPYba, 88 8b,dPPYba,   ,adPPYb,d8 88  ,adPPYba,    8b,dPPYba,   ,adPPYba, 8b,dPPYba,   
I8[    "" 88 88P'   `"8a a8"    `Y88 88 a8P_____88    88P'    "8a a8P_____88 88P'   `"8a  
 `"Y8ba,  88 88       88 8b       88 88 8PP"""""""    88       d8 8PP""""""" 88       88  
aa    ]8I 88 88       88 "8a,   ,d88 88 "8b,   ,aa    88b,   ,a8" "8b,   ,aa 88       88  
`"YbbdP"' 88 88       88  `"YbbdP"Y8 88  `"Ybbd8"'    88`YbbdP"'   `"Ybbd8"' 88       88
                          aa,    ,88                  88
                           "Y8bbdP"                   88

*/

//simulates the motion of a single pendulum
//using multiple finite difference methods
void single_pendulum()
{
	char processName[64] = "Single Pendulum";

	ofstream single_pen("data/single_pen_combined.csv");

	if(outputPositions)
	{
		single_pen << ",euler_theta,euler_w" ;
		single_pen << ",leapfrog_theta,leapfrog_w" ;
		single_pen << ",rk4_theta,rk4_w" ;
	}

	if(outputEnergies)
	{
		single_pen << ",euler_T,euler_U, euler_E" ;
		single_pen << ",leapfrog_T,leapfrog_U, leapfrog_E" ;
		single_pen << ",rk4_T,rk4_U, rk4_E" ;
	}
	//end column headers
	single_pen << "\n";

	ofstream euler_single_pen("data/euler_single_pen.csv");
	ofstream leapfrog_single_pen("data/leapfrog_single_pen.csv");
	ofstream rk4_single_pen("data/rk4_single_pen.csv");

	initFile(euler_single_pen, "euler");
	initFile(leapfrog_single_pen, "leapfrog");
	initFile(rk4_single_pen, "rk4");

	//sets the number of decimal places time is outputted to
	//see: http://www.cplusplus.com/reference/ios/scientific/
	std::cout.precision(1);


/*
          88                           88                                      
          ""                           88                                      
                                       88                                      
,adPPYba, 88 88,dPYba,,adPYba,         88  ,adPPYba,   ,adPPYba,  8b,dPPYba,   
I8[    "" 88 88P'   "88"    "8a        88 a8"     "8a a8"     "8a 88P'    "8a  
 `"Y8ba,  88 88      88      88        88 8b       d8 8b       d8 88       d8  
aa    ]8I 88 88      88      88 888    88 "8a,   ,a8" "8a,   ,a8" 88b,   ,a8"  
`"YbbdP"' 88 88      88      88 888    88  `"YbbdP"'   `"YbbdP"'  88`YbbdP"'   
                                                                  88           
                                                                  88 
*/

	for (int i = 0; i < numberOfSteps; i++)
	{
		updateEnergies(i);

		/************** OUTPUT TO FILE **************/

		//output time
		t = h*i;
		euler_single_pen << std::fixed << t ;
		leapfrog_single_pen << std::fixed << t ;
		rk4_single_pen << std::fixed << t ;
	

		if(outputPositions)
		{
			outputPositionToFile(euler_single_pen, euler_theta[i], euler_w[i]);
			outputPositionToFile(leapfrog_single_pen, leapfrog_theta[i], leapfrog_w[i]);
			outputPositionToFile(rk4_single_pen, rk4_theta[i], rk4_w[i]);
		}

		if(outputEnergies)
		{
			outputEnergyToFile(euler_single_pen, euler_T[i], euler_U[i], euler_E[i]);
			outputEnergyToFile(leapfrog_single_pen, leapfrog_T[i], leapfrog_U[i], leapfrog_E[i]);
			outputEnergyToFile(rk4_single_pen, rk4_T[i], rk4_U[i], rk4_E[i]);
		}

		euler_single_pen << "\n" ;
		leapfrog_single_pen << "\n" ;
		rk4_single_pen << "\n" ;
		
		/****************************************************************/
		//output time
		single_pen << std::fixed << t ;

		if(outputPositions)
		{
			single_pen << std::scientific << "," << euler_theta[i] << "," << euler_w[i] ;
			single_pen << std::scientific << "," << leapfrog_theta[i] << "," << leapfrog_w[i] ;
			single_pen << std::scientific << "," << rk4_theta[i] << "," << rk4_w[i] ;
		}

		if(outputEnergies)
		{
			single_pen << std::scientific << "," << euler_T[i] << "," << euler_U[i] << "," << euler_E[i];
			single_pen << std::scientific << "," << leapfrog_T[i] << "," << leapfrog_U[i] << "," << leapfrog_E[i];
			single_pen << std::scientific << "," << rk4_T[i] << "," << rk4_U[i] << "," << rk4_E[i];
		}

		//end output for this iteration
		single_pen << "\n";
		/****************************************************************/
		
		updateEuler(i);
		updateLeapfrog(i);
		updateRK4(i);

		//update progress meter in terminal
		updateProgress(i, processName);

	}
} //end single_pendulum()

void updateEnergies(int i)
{
	//euler
	euler_T[i] = calculateKineticEnergy(m, l, euler_theta[i], euler_w[i]);
	euler_U[i] = calculatePotentialEnergy(m, l, euler_theta[i], euler_w[i]);
	euler_E[i] = euler_T[i] + euler_U[i];

	//leapfrog
	leapfrog_T[i] = calculateKineticEnergy(m, l, leapfrog_theta[i], leapfrog_w[i]);
	leapfrog_U[i] = calculatePotentialEnergy(m, l, leapfrog_theta[i], leapfrog_w[i]);
	leapfrog_E[i] = leapfrog_T[i] + leapfrog_U[i];

	//RK4
	rk4_T[i] = calculateKineticEnergy(m, l, rk4_theta[i], rk4_w[i]);
	rk4_U[i] = calculatePotentialEnergy(m, l, rk4_theta[i], rk4_w[i]);
	rk4_E[i] = rk4_T[i] + rk4_U[i];
}

void outputPositionToFile(ostream& file, double theta, double w)
{
	file << std::scientific << "," << theta << "," << w ;
}

void outputEnergyToFile(ostream& file, double T, double U, double E)
{
	file << std::scientific << "," << T << "," << U << "," << E;
}

void updateEuler(int i)
{
	//update euler_theta
	euler_theta[i+1] =  ( h*euler_w[i] ) + euler_theta[i];
	//update euler_w
	euler_w[i+1] = -h*euler_theta[i] + ( 1 - h*beta )*euler_w[i];
}

void updateLeapfrog(int i)
{
	if(i==0) //we need to use the euler method to give us some initial values
	{
		//update leapfrog_theta
		leapfrog_theta[i+1] =  euler_theta[i+1];
		//update leapfrog_w
		leapfrog_w[i+1] = euler_w[i+1];
	}
	else
	{
		//update leapfrog_theta
		leapfrog_theta[i+1] =  leapfrog_theta[i-1] + leapfrog_w[i]*2.0*h;
		//update leapfrog_w
		leapfrog_w[i+1] = leapfrog_w[i-1] - 2.0*h*( leapfrog_theta[i] + ( beta*leapfrog_w[i] ) );
	}
}

void updateRK4(int i)
{
	k_1 = h * ( rk4_w[i] );
	k_2 = h * ( rk4_w[i] + 0.5*k_1 );
	k_3 = h * ( rk4_w[i] + 0.5*k_2 );
	k_4 = h * ( rk4_w[i] + k_3);

	rk4_theta[i+1] = rk4_theta[i] + (1.0/6.0)*(k_1 + 2*k_2 + 2*k_3 + k_4);

	k_1 = -1*h * ( rk4_theta[i] + beta*rk4_w[i] );
	k_2 = -1*h * ( rk4_theta[i]+0.5*h + beta*(rk4_w[i] + 0.5*k_1) );
	k_3 = -1*h * ( rk4_theta[i]+0.5*h + beta*(rk4_w[i] + 0.5*k_2) );
	k_4 = -1*h * ( rk4_theta[i]+h + beta*(rk4_w[i] + k_3) );

	rk4_w[i+1] = rk4_w[i] + (1.0/6.0)*(k_1 + 2*k_2 + 2*k_3 + k_4);
}

void stabilityTest() {
	calculateAnalyticalSolution();
	
	double anal_E_initial	= calculateKineticEnergy(m, l, anal_theta[0], anal_w[0]);
	double anal_E_final		= calculateKineticEnergy(m, l, anal_theta[numberOfSteps-1], anal_w[numberOfSteps-1]);

	string euler_stable			= isStable(euler_E[0], euler_E[numberOfSteps-1]);
	string leapfrog_stable	= isStable(leapfrog_E[0], leapfrog_E[numberOfSteps-1]);
	string rk4_stable				= isStable(rk4_E[0], rk4_E[numberOfSteps-1]);

	ofstream stabilityTestLog("data/stability_tests.csv");
	stabilityTestLog << "Method,Initial Energy,Final Energy,Stable?\n" ;

	stabilityTestLog << "Euler," 		<< euler_E[0]			<< "," << euler_E[numberOfSteps-1]		<< "," << euler_stable		<< "\n";
	stabilityTestLog << "Leapfrog,"	<< leapfrog_E[0]	<< "," << leapfrog_E[numberOfSteps-1]	<< "," << leapfrog_stable	<< "\n";
	stabilityTestLog << "RK4,"			<< rk4_E[0]				<< "," << rk4_E[numberOfSteps-1]			<< "," << rk4_stable			<< "\n";
}

string isStable(double E_i, double E_f) {
	if(E_f > E_i)
	{
		return "stable";
	}
	else
	{
		return "unstable";
	}
}

void eulerStabilityTest()
{
	ofstream stabilityTestLog("data/euler_stability_test.csv");
	calculateAnalyticalSolution();
	stabilityTestLog << "time,euler_theta,anal_theta,g,stable?\n" ;
	double g_stab = 0.0;
	for (int i = 0; i < numberOfSteps; ++i)
	{
		err_theta[i] = euler_theta[i] - anal_theta[i];
		if(i != 0)
		{
			g_stab = abs(err_theta[i]/err_theta[i-1]);
			g_stab = floorf(g_stab * 10000 + 0.5) / 10000; //crude rounding method
		}
		if(g_stab<1.0)
		{
			stabilityTestLog << i*h << "," << euler_theta[i] << "," << anal_theta[i] << "," << g_stab << "," << "HELL YES!\n";
		}
		else
		{
			stabilityTestLog << i*h << "," << euler_theta[i] << "," << anal_theta[i] << "," << g_stab << "," << "NOPE\n";
		}
	}
}

void calculateAnalyticalSolution()
{
	for (int i = 0; i < numberOfSteps; ++i)
	{
		anal_theta[i] = initial_theta*cos(i*h);
	}
}

/*
                                                                       
88                     88                                              
88                     88                                              
88                     88                                              
88,dPPYba,   ,adPPYba, 88 8b,dPPYba,   ,adPPYba, 8b,dPPYba, ,adPPYba,  
88P'    "8a a8P_____88 88 88P'    "8a a8P_____88 88P'   "Y8 I8[    ""  
88       88 8PP""""""" 88 88       d8 8PP""""""" 88          `"Y8ba,   
88       88 "8b,   ,aa 88 88b,   ,a8" "8b,   ,aa 88         aa    ]8I  
88       88  `"Ybbd8"' 88 88`YbbdP"'   `"Ybbd8"' 88         `"YbbdP"'  
                          88                                           
                          88 
*/

//calculate kinetic energy
double calculateKineticEnergy(double m, double l, double theta, double w)
{
	return 0.5*m*pow(l,2.0)*pow(w,2.0);
}

//calculate potential energy
double calculatePotentialEnergy(double m, double l, double theta, double w)
{
	return m*g*l*(1 - cos(theta));
}

//updates the progress display on commandline
void updateProgress(int i, char *name)
{
		cout << std::fixed << "\r"<< "Currently running " << name << " - " << (float(i)/numberOfSteps)*100 << "\%" << flush;
}

//displays done message to terminal
void done()
{
	std::cout << "\r"<< "100\% - all done. Simulated time elapsed: " << simulatedTime << "s\n" << std::flush;
}