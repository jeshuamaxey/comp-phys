#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <string>

#define simulatedTime 60.0												//simulated time (seconds)
#define h_min 0.01
#define h_max 0.5
#define h_step 0.02
#define numberOfSteps int(simulatedTime / h_min)	//used for sizing arrays
#define g 9.81																		//acceleration due to gravity
#define pi atan(1.0)															//mutha fuckin pi man

#define outputPositions true											//determine what's outputted by the program to file
#define outputEnergies false											//determine what's outputted by the program to file

using namespace std;

//double pendulum functions
void double_pendulum(double, double, double);
void updateRK4(double, double, double, int);
string makeFileName(double h, double G, double R);

//helper functions
void updateProgress(int);
void updateProgress(int, char[]);
void done();

int main()
{
	char processName[64] = "Double Pendulum";
	
	double h;
	double G = 0;																			//matrix constant - gamma/(m* sqrt( g*l ))
	double R = 100;																		//matrix constant - M/m

	double G_min = 0;
	double G_max = 1;
	double G_step = 1;

	double R_min = 1;
	double R_max = 101;
	double R_step = 10;

	int h_range = int( (h_max - h_min)/h_step );
	int G_range = int( (G_max - G_min)/G_step );
	int R_range = int( (R_max - R_min)/R_step );

	//yo dawg I heard you like for loops
	for (int i = 0; i < h_range; ++i)
	{
		updateProgress(i, processName);
		h = h_min + i*h_step;
		for (int i = 0; i < G_range; ++i)
		{
			G = G_min + i*G_step;
			for (int i = 0; i < R_range; ++i)
			{
				R = R_min + i*R_step;
				//make the call
				double_pendulum(h, G, R);
			}
		}
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
void double_pendulum(double h, double G, double R)
{

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

	string fileName = makeFileName(h, G, R);
	//open up file and set column headings
	ofstream double_pen(fileName);
	//first column is always time
	double_pen << "time" ;
	//rk4
	double_pen << ",rk4_theta,rk4_psi,rk4_w,rk4_v" ;
	//end column headers
	double_pen << "\n";

	for (int i = 0; i < numberOfSteps; i++)
	{
		/************** OUTPUT **************/
		//output time
		double_pen << std::fixed << h*i ;
		//output rk4
		double_pen << std::scientific << "," << rk4_theta[i] << "," << rk4_psi[i] << "," << rk4_w[i] << "," << rk4_v[i] ;
		//end output for this iteration
		double_pen << "\n";

		updateRK4(h, R, G, i);

	}
} //end double_pendulum

void updateRK4(double h, double R, double G, int i)
{
	/************** NEW RK4 METHOD **************/
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
	name << "data/dp/double_pen_h="<< h << "_G=" << G << "_R=" << R << ".csv";
	return name.str();
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