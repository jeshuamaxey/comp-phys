#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <string>

#define simulatedTime 50.0											//simulated time (seconds)
#define h_min 0.01
#define h_max 4.0
#define h_step 0.005
#define numberOfSteps int(simulatedTime / h_min)	//used for sizing arrays
#define g 9.81																	//acceleration due to gravity
#define pi atan(1.0)														//mutha fuckin pi man

#define outputPositions true										//determine what's outputted by the program to file
#define outputEnergies true											//determine what's outputted by the program to file
#define runStabilityTest true

using namespace std;

//returns the number of loop iterations it ran for a particular simulation
int single_pendulum(double, double);
void updateStabilityTestFile(ostream&, double, double, int);
void setInitialValues(double, double);

string makeFileName(double, double, string);
string makeStabTestFileName(string);
void initFile(ostream&, string);
string makePositionHeading(string);
string makeEnergyHeading(string);

void outputPositionToFile(ostream&, double, double);
void outputEnergyToFile(ostream&, double);

void updateEuler(int, double, double);
void updateLeapfrog(int, double, double);
void updateRK4(int, double, double);
void updateEnergies(int, double, double);

//helper functions
void updateProgress(double, char[]);
void done();

/********** GLOBAL VARIABLES COS FUCK SCOPE ***********/

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
double k_1[2], k_2[2], k_3[2], k_4[2];

//analytical solution & stabiliy testvars
double anal_theta[numberOfSteps];
double anal_w[numberOfSteps];
double err_theta[numberOfSteps];

/***************************************************************/


int main()
{
	char processName[64] = "Single Pendulum";

	double h, damping_constant;

	double damping_constant_min = 0.2;
	double damping_constant_max = 0.2;
	double damping_constant_step = 0.2;

	int h_range = int( (h_max - h_min)/h_step );
	int damping_constant_range = int( (damping_constant_max - damping_constant_min)/damping_constant_step );

	//create file to output stability test data to
	ofstream energyStabTestFile("data/sp/energy_stab_test/data.csv");
	//set column headings
	energyStabTestFile << "h,Euler E_final/E_initial,Leapfrog E_final/E_initial,RK4 E_final/E_initial\n";

	for (int i = 0; i < h_range; i++)
	{
		updateProgress(float(i)/h_range, processName);
		h = h_min + i*h_step;
		for (int j = 0; j <= damping_constant_range; j++)
		{
			damping_constant = damping_constant_min + j*damping_constant_step;
			//make the call
			int maxIndex = single_pendulum(h, damping_constant);
			
			if(runStabilityTest)
			{
				energyStabTestFile << h ;
				updateStabilityTestFile(energyStabTestFile, euler_E[1], euler_E[maxIndex], maxIndex);
				updateStabilityTestFile(energyStabTestFile, leapfrog_E[1], leapfrog_E[maxIndex], maxIndex);
				updateStabilityTestFile(energyStabTestFile, rk4_E[1], rk4_E[maxIndex], maxIndex);
				energyStabTestFile << "\n";
			}
		}
	}
	done();
	return 0;
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
int single_pendulum(double h, double damping_constant)
{
	double t = 0;																						//time
	double l = 9.81;																				//length of pendulem in metres
	double m = 1.0;																					//mass of pendulum in kg
	double beta = damping_constant/(m* sqrt( g*l ));				//matrix constant

	//initial conditions
	double initial_theta = 0.01;															//angle from vert (starting angle)
	double initial_w = 0;																		//rate of change of angle from vert (starting at rest)

	setInitialValues(initial_theta, initial_w);
	updateEnergies(0, m, l);

	std::cout.precision(1);																	//sets the number of decimal places time is outputted to
																													//see: http://www.cplusplus.com/reference/ios/scientific/

	//create custom filenames
	string eulerFileName		= makeFileName(h, damping_constant, "euler");
	string leapfrogFileName	= makeFileName(h, damping_constant, "leapfrog");
	string rk4FileName			= makeFileName(h, damping_constant, "rk4");
	//open up files
	ofstream eulerFile(eulerFileName);
	ofstream leapfrogFile(leapfrogFileName);
	ofstream rk4File(rk4FileName);
	//set column headings
	initFile(eulerFile, "euler");
	initFile(leapfrogFile, "leapfrog");
	initFile(rk4File, "rk4");

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

  //set loop limit
  int loopLimit = int(simulatedTime / h);
	for (int i = 0; i < loopLimit; i++)
	{
		

		/************** OUTPUT TO FILE **************/

		//output time
		t = h*i;
		eulerFile << std::fixed << t ;
		leapfrogFile << std::fixed << t ;
		rk4File << std::fixed << t ;
	

		if(outputPositions)
		{
			outputPositionToFile(eulerFile, euler_theta[i], euler_w[i]);
			outputPositionToFile(leapfrogFile, leapfrog_theta[i], leapfrog_w[i]);
			outputPositionToFile(rk4File, rk4_theta[i], rk4_w[i]);
		}

		if(outputEnergies)
		{
			outputEnergyToFile(eulerFile, euler_E[i]);
			outputEnergyToFile(leapfrogFile, leapfrog_E[i]);
			outputEnergyToFile(rk4File, rk4_E[i]);
		}

		eulerFile << "\n" ;
		leapfrogFile << "\n" ;
		rk4File << "\n" ;
		
		updateEuler(i, h, beta);
		updateLeapfrog(i, h, beta);
		updateRK4(i, h, beta);
		updateEnergies(i+1, m, l);

	}
	return loopLimit-1;
} //end single_pendulum()

void setInitialValues(double theta, double w)
{
	//set initial values
	euler_theta[0] = theta;
	euler_w[0] = w;

	leapfrog_theta[0] = theta;
	leapfrog_w[0] = w;

	rk4_theta[0] = theta;
	rk4_w[0] = w;
}

string makeFileName(double h, double gamma, string processName)
{
	stringstream name;
	name << "data/sp/single_pen_" << processName << "_h=" << h << "_gamma=" << gamma << ".csv";
	return name.str();
}

string makeStabTestFileName(string processName)
{
	stringstream name;
	name << "data/sp/energy_stab_tests/" << processName << ".csv";
	return name.str();
}

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
	string << "," << processName << "_E" ;
	return string.str();
}

void updateEnergies(int i, double m, double l)
{
	//euler
	euler_E[i] = pow(euler_theta[i], 2.0) + pow(euler_w[i], 2.0);
	//leapfrog
	leapfrog_E[i] = pow(leapfrog_theta[i], 2.0) + pow(leapfrog_w[i], 2.0);
	//RK4
	rk4_E[i] = pow(rk4_theta[i], 2.0) + pow(rk4_w[i], 2.0);
}

void outputPositionToFile(ostream& file, double theta, double w)
{
	file << std::scientific << "," << theta << "," << w ;
}

void outputEnergyToFile(ostream& file, double E)
{
	file << std::scientific << "," << E;
}

void updateEuler(int i, double h, double beta)
{
	//update euler_theta
	euler_theta[i+1] =  ( h*euler_w[i] ) + euler_theta[i];
	//update euler_w
	euler_w[i+1] = -h*euler_theta[i] + ( 1 - h*beta )*euler_w[i];
}

void updateLeapfrog(int i, double h, double beta)
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

void updateRK4(int i, double h, double beta)
{

	k_1[0] = h * ( rk4_w[i] );
	k_1[1] = -1*h * ( rk4_theta[i] + beta*rk4_w[i] );

	k_2[0] = h * ( rk4_w[i] + 0.5*k_1[1] );
	k_2[1] = -1*h * ( rk4_theta[i]+0.5*k_1[0] + beta*(rk4_w[i] + 0.5*k_1[1]) );

	k_3[0] = h * ( rk4_w[i] + 0.5*k_2[1] );
	k_3[1] = -1*h * ( rk4_theta[i]+0.5*k_2[0] + beta*(rk4_w[i] + 0.5*k_2[1]) );

	k_4[0] = h * ( rk4_w[i] + k_3[1]);
	k_4[1] = -1*h * ( rk4_theta[i]+k_3[0] + beta*(rk4_w[i] + k_3[1]) );


	rk4_theta[i+1] = rk4_theta[i] + (1.0/6.0)*(k_1[0] + 2*k_2[0] + 2*k_3[0] + k_4[0]);
	rk4_w[i+1] = rk4_w[i] + (1.0/6.0)*(k_1[1] + 2*k_2[1] + 2*k_3[1] + k_4[1]);
}

void updateStabilityTestFile(ostream& file, double E_i, double E_f, int i)
{
	double ratio = E_f / E_i;
	//ratio != ratio checks to see if ratio == nan
	if(ratio > 4.0 || ratio != ratio)
	{
		file << "," <<  4.0;
	}
	else
	{
		file << "," <<  ratio;
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

//updates the progress display on commandline
void updateProgress(double progress, char *name)
{
		cout << std::fixed << "\r" << "Currently running " << name << " - " << int(progress*100) << "\%" << flush;
}

//displays done message to terminal
void done()
{
	std::cout << "\r"<< "100\% - all done. Simulated time elapsed: " << simulatedTime << "s\n" << std::flush;
}