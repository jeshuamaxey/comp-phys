#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <string>

#define simulatedTime 100.0												//simulated time (seconds)
#define h_min 0.01
#define h_max 0.02
#define h_step 0.01
#define numberOfSteps int(simulatedTime / h_min)	//used for sizing arrays
#define g 9.81																		//acceleration due to gravity
#define pi atan(1.0)															//mutha fuckin pi man

#define runAutoTests true													//should be set to true most of the time
#define runSpecificTests false

//recommend seeting only one of these to true at a time
//they all output to the same filestream and while column
//headings should be correct, there's no guaurantee
#define outputPositions true											//determine what's outputted by the program to file
#define outputEnergies false											//determine what's outputted by the program to file
#define outputIndividualPendulumEnergies false
#define runStabTest false													//this takes a long time - set to true at your peril (was worth it for the graph though)

//for debug purposes or if you don't want to overwrite data
#define useTestDir false

using namespace std;

//double pendulum functions
void double_pendulum(double, double, double, ostream&);
void updateRK4(double, double, double, int);
double calculateTotalEnergy(double, double, double, double, double, double);
double calculatePotentialEnergy(double, double, double, double, double, double);
double calculateKineticEnergy(double, double, double, double, double, double);
double calculateEnergyHigher(double theta, double psi, double w, double v, double R, double G );
double calculateEnergyLower(double theta, double psi, double w, double v, double R, double G );
string makeFileName(double h, double G, double R);

//helper functions
void updateProgress(double, char[]);
void done();

//global debug file stream object
ofstream debug("data/__debug_log.txt");
//global stab test file stream object (couldn't get it to work if it wasn't global...)
ofstream stab_test_output("data/dp/_stab_test_results/stab_test.csv");
//used in stability test - global cos I'm in a hurry here ok?
double h_prev = h_min;

int main()
{
	stab_test_output << "h,energy ratio R=0.001,energy ratio R=0.01,energy ratio R=0.1,energy ratio R=1,energy ratio R=10,energy ratio R=100,energy ratio R=1000\n" ;

	char processName[64] = "Double Pendulum";
	
	double h;

	if(runAutoTests)
	{
		double G, R;
		double G_min = 0;
		double G_max = 1;
		double G_step = 1;

		double R_min = 0.001;
		double R_max = 1000;
		double R_step = 10;
		
		double R_scale = 10.0;

		int h_range = int( (h_max - h_min)/h_step );
		int G_range = int( (G_max - G_min)/G_step );

		//used when R is increased linearly
		//int R_range = int( (R_max - R_min)/R_step );
		
		//used when R is increased by scale factor
		//int R_range = int( (R_max/R_min) / R_scale );
		int R_range = 6;
		
		debug << "h_range = " << h_range << "\n";

		//yo dawg I heard you like for loops
		for (int i = 0; i < h_range; ++i)
		{
			updateProgress(float(i)/h_range, processName);
			h = h_min + i*h_step;
			for (int i = 0; i <= G_range; ++i)
			{
				G = G_min + i*G_step;
				for (int j = 0; j <= R_range; ++j)
				{
					R = R_min * pow(R_scale, j);

					//the stab test file requires a different column format
					//and creation before each simulation begins
					stringstream stab_test_file_name;
					stab_test_file_name << "data/dp/_stab_test_results/R=" << R << "_G=" << G << ".csv" ;
					ofstream double_pen_stab_test(stab_test_file_name.str());
					double_pen_stab_test << "h,E ratio\n" ;

					//debug << "h = " << h << "\n";
					//make the call
					double_pendulum(h, G, R, double_pen_stab_test);
				} //end R for loop
			} //end G for loop
		} //end h for loop
	}

	if(runSpecificTests)
	{
		int G_arr[2] = {0, 1};
		double R_arr[3] = {0.01, 0.1, 100.0};

		int h_range = int( (h_max - h_min)/h_step );

		//yo dawg I heard you like for loops
		for (int i = 0; i <= h_range; ++i)
		{
			updateProgress(float(i)/h_range, processName);
			h = h_min + i*h_step;
			for (int j = 0; j < 2; ++j)
			{
				for (int k = 0; k < 3; ++k)
				{
					//dummy file, not used
					ofstream dummy("DELETE_ME");
					//make the call
					double_pendulum(h, G_arr[j], R_arr[k], dummy);
				} //end R for loop
			} //end G for loop
		} //end h for loop

	}
	done();
	return 0;
}

/*
                                                      
         88                         88          88             
         88                         88          88             
         88                         88          88             
 ,adPPYb,88  ,adPPYba,  88       88 88,dPPYba,  88  ,adPPYba,  
a8"    `Y88 a8"     "8a 88       88 88P'    "8a 88 a8P_____88  
8b       88 8b       d8 88       88 88       d8 88 8PP"""""""  
"8a,   ,d88 "8a,   ,a8" "8a,   ,a88 88b,   ,a8" 88 "8b,   ,aa  
 `"8bbdP"Y8  `"YbbdP"'   `"YbbdP'Y8 8Y"Ybbd8"'  88  `"Ybbd8"'  
                                                               
                                                               
                                    
8b,dPPYba,   ,adPPYba, 8b,dPPYba,   
88P'    "8a a8P_____88 88P'   `"8a  
88       d8 8PP""""""" 88       88  
88b,   ,a8" "8b,   ,aa 88       88  
88`YbbdP"'   `"Ybbd8"' 88       88  
88                                  
88   

*/

/********************************* GLOBAL VARIABLES ********************************/

double rk4_theta [numberOfSteps];						//angular displacement of 1st pendulum
double rk4_psi [numberOfSteps];							//angular displacement of 2nd pendulum
double rk4_w [numberOfSteps];								//angular acceleration of 1st pendulum
double rk4_v [numberOfSteps];								//angular acceleration of 2nd pendulum

double k1[4], k2[4], k3[4], k4[4];					//RK4 values

/************************************************************************************/

//simulates the motion of a double pendulum
//using the RK4 finite difference method
void double_pendulum(double h, double G, double R, ostream& stab_test_file)
{
	//debug << "h = " << h << "\n";

	double initial_theta = 0.1;									//angle from vert for pendulum 1 (starting angle)
	double initial_psi = 0.0;										//angle from vert for pendulum 2 (starting angle)
	double initial_w = 0.0;											//rate of change of angle from vert for pendulum 1 (starting at rest)
	double initial_v = 0.0;											//rate of change of angle from vert for pendulum 2 (starting at rest)

	std::cout.precision(1);											//sets the number of decimal places time is outputted to
																							//see: http://www.cplusplus.com/reference/ios/scientific/

	//set initial values
	rk4_theta[0] = initial_theta;
	rk4_psi[0] = initial_psi;
	rk4_w[0] = initial_w;
	rk4_v[0] = initial_v;

	double initial_energy = calculateTotalEnergy(rk4_theta[0], rk4_psi[0], rk4_w[0], rk4_v[0], R, G);

	//same file used to output all data
	string fileName = makeFileName(h, G, R);
	//open up file and set column headings
	ofstream double_pen(fileName);

	//first column is always time
	double_pen << "time" ;

	if(outputPositions)
	{
		double_pen << ",rk4_theta,rk4_psi,rk4_w,rk4_v" ;
	}
	if(outputEnergies)
	{
		double_pen << ",U,T,E" ;
	}
	if(outputIndividualPendulumEnergies)
	{
		double_pen << ",E_higher,E_lower" ;
	}
	//end column headers
	double_pen << "\n";

	//find number of iterations required for this loop
	//not used for simulation loop - used for stab test
	int maxIndex = int( simulatedTime / h );

	//debug << "maxIndex = " << maxIndex << "\n";

	for (int i = 0; i < numberOfSteps; i++)
	{
		/************** OUTPUT **************/
		//output time
		double_pen << std::fixed << h*i ;
		
		if(outputPositions)
		{
			double_pen << std::scientific << "," << rk4_theta[i] << "," << rk4_psi[i] << "," << rk4_w[i] << "," << rk4_v[i] ;
		}
		if(outputEnergies)
		{
			double U = calculatePotentialEnergy(rk4_theta[i], rk4_psi[i], rk4_w[i], rk4_v[i], R, G );
			double T = calculateKineticEnergy(rk4_theta[i], rk4_psi[i], rk4_w[i], rk4_v[i], R, G );
			double_pen << std::scientific << "," << U
																		<< "," << T
																		<< "," << U+T ;
		}
		if(outputIndividualPendulumEnergies)
		{
			double_pen << std::scientific << "," << calculateEnergyHigher(rk4_theta[i], rk4_psi[i], rk4_w[i], rk4_v[i], R, G)
																		<< "," << calculateEnergyLower(rk4_theta[i], rk4_psi[i], rk4_w[i], rk4_v[i], R, G) ;
		}
		if(runStabTest && (i == int(maxIndex-1)))
		{
			//calculate E_final / E_initial
			double ratio = calculateTotalEnergy(rk4_theta[i], rk4_psi[i], rk4_w[i], rk4_v[i], R, G) / calculateTotalEnergy(rk4_theta[0], rk4_psi[0], rk4_w[0], rk4_v[0], R, G );
			
			//if ratio is greater than 4 or not a number, set it to 4 to aid graph plotting
			if(ratio > 1.5 || ratio != ratio)
			{
				ratio = 1.5;
			}
			//debug << "i=" << i << ", maxIndex=" << maxIndex << ", h=" << h << ", ratio=" << ratio << "\n";
			//debug << h << "," << G << "," << R << "," << ratio << "\n";
			if(h != h_prev)
			{
				stab_test_output << "\n" << h;
				h_prev = h;
			}
			stab_test_output << "," << ratio;
		}
		//end output for this iteration
		double_pen << "\n";

		updateRK4(h, R, G, i);
	}
} //end double_pendulum

void updateRK4(double h, double R, double G, int i)
{
	/************** NEW RK4 METHOD **************/
	//time permitting I might refactor this to
	//loop through an update matrix, but it works
	//in its current form.
	k1[0] = h * rk4_w[i];
	k1[1] = h * rk4_v[i];
	k1[2] = h * ( -(R+1.0)*rk4_theta[i] + R*rk4_psi[i] - G*rk4_w[i] );
	k1[3] = h * ( (R+1.0)*rk4_theta[i] -(R+1.0)*rk4_psi[i] + G*(1.0-(1.0/R))*rk4_w[i] - (G/R)*rk4_v[i] );

	k2[0] = h * ( rk4_w[i] 		+ 0.5*k1[2] );
	k2[1] = h * ( rk4_v[i]	+ 0.5*k1[3] );
	k2[2] = h * ( -(R+1.0)*(rk4_theta[i] + 0.5*k1[0]) + R*(rk4_psi[i] + 0.5*k1[1]) - G*(rk4_w[i] + 0.5*k1[2]) );
	k2[3] = h * ( (R+1.0)*(rk4_theta[i] + 0.5*k1[0]) -(R+1.0)*(rk4_psi[i] + 0.5*k1[1]) + G*(1.0-(1.0/R))*(rk4_w[i] + 0.5*k1[2]) - (G/R)*(rk4_v[i] + 0.5*k1[3]) );

	k3[0] = h * ( rk4_w[i] 		+ 0.5*k2[2] );
	k3[1] = h * ( rk4_v[i]	+ 0.5*k2[3] );
	k3[2] = h * ( -(R+1.0)*(rk4_theta[i] + 0.5*k2[0]) + R*(rk4_psi[i] + 0.5*k2[1]) - G*(rk4_w[i] + 0.5*k2[2]) );
	k3[3] = h * ( (R+1.0)*(rk4_theta[i] + 0.5*k2[0]) -(R+1.0)*(rk4_psi[i] + 0.5*k2[1]) + G*(1.0-(1.0/R))*(rk4_w[i] + 0.5*k2[2]) - (G/R)*(rk4_v[i] + 0.5*k2[3]) );

	k4[0] = h * ( rk4_w[i] 		+ k3[2] );
	k4[1] = h * ( rk4_v[i]	+ k3[3] );
	k4[2] = h * ( -(R+1.0)*(rk4_theta[i] + k3[0]) + R*(rk4_psi[i] + k3[1]) - G*(rk4_w[i] + k3[2]) );
	k4[3] = h * ( (R+1.0)*(rk4_theta[i] + k3[0]) -(R+1.0)*(rk4_psi[i] + k3[1]) + G*(1.0-(1.0/R))*(rk4_w[i] + k3[2]) - (G/R)*(rk4_v[i] + k3[3]) );

	rk4_theta[i+1]	= rk4_theta[i]	+ (1.0/6.0)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
	rk4_psi[i+1]		= rk4_psi[i]		+ (1.0/6.0)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
	rk4_w[i+1]			= rk4_w[i]			+ (1.0/6.0)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
	rk4_v[i+1]			= rk4_v[i]			+ (1.0/6.0)*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
}

//returns total energy of system
double calculateTotalEnergy(double theta, double psi, double w, double v, double R, double G )
{
	return calculatePotentialEnergy(theta, psi, w, v, R, G )
					+ calculateKineticEnergy(theta, psi, w, v, R, G );
}

//returns potential energy of system
double calculatePotentialEnergy(double theta, double psi, double w, double v, double R, double G )
{
	double M = R;
	double m = 1.0;
	double l = 1.0;

	//
	double adjustment = sqrt(g/l);

	w *= adjustment;
	v *= adjustment;

	return 0.5*g*( (m+M)*l*theta*theta + M*l*psi*psi );
}

//returns kinetic energy of system
double calculateKineticEnergy(double theta, double psi, double w, double v, double R, double G )
{
	double M = R;
	double m = 1.0;
	double l = 1.0;

	double adjustment = sqrt(g/l);

	w *= adjustment;
	v *= adjustment;

	return 0.5*pow(l, 2.0)*( m*pow(w, 2.0) + M*( pow(w, 2.0) + pow(v, 2.0) +2*w*v ) );
}

//returns energy of the higher pendulum bob
double calculateEnergyHigher(double theta, double psi, double w, double v, double R, double G )
{
	double M = R;
	double m = 1.0;
	double l = 1.0;

	//
	double adjustment = sqrt(g/l);

	w *= adjustment;
	v *= adjustment;

	return 0.5*g*(m+M)*l*theta*theta + 0.5*pow(l, 2.0)*m*pow(w, 2.0);
}

//returns energy of the lower pendulum bob
double calculateEnergyLower(double theta, double psi, double w, double v, double R, double G )
{
	double M = R;
	double m = 1.0;
	double l = 1.0;

	double adjustment = sqrt(g/l);

	w *= adjustment;
	v *= adjustment;

	return 0.5*pow(l, 2.0)*( M*( pow(w, 2.0) + pow(v, 2.0) +2*w*v ) ) + 0.5*g*M*l*psi*psi;
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

string makeFileName(double h, double G, double R)
{
	stringstream name;
	if(useTestDir)
	{
		name << "data/dp_test/double_pen_h="<< h << "_G=" << G << "_R=" << R << ".csv";
	}
	else
	{
		name << "data/dp/double_pen_h="<< h << "_G=" << G << "_R=" << R << ".csv";
	}
	return name.str();
}

//updates the progress display on commandline
void updateProgress(double progress, char *name)
{
		std::cout << "\r"<< "Currently running " << name << " - " << int(progress*100) << "%" << std::flush;
}

//displays done message to terminal
void done()
{
	std::cout << std::fixed << "\r"<< "100\% - all done. Simulated time elapsed: " << simulatedTime << "s\n" << std::flush;
}