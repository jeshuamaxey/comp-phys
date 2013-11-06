#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#define simulatedTime 6.0											//simulated time (seconds)
#define h 0.01																//step size
#define numberOfSteps  int(simulatedTime / h)	//used for sizing arrays
#define g 9.81																//acceleration due to gravity
#define pi atan(1.0)													//mutha fuckin pi man

#define outputPositions true									//determine what's outputted by the program to file
#define outputEnergies false									//determine what's outputted by the program to file

using namespace std;


//single pendulum functions
void single_pendulum();
double calculateKineticEnergy(double m, double l, double theta, double w);
double calculatePotentialEnergy(double m, double l, double theta, double w);

//double pendulum functions
void double_pendulum();

//helper functions
void updateProgress(int);
void updateProgress(int, char[]);
void done();

int main()
{
	single_pendulum();
	double_pendulum();
	done();
	return 0;
}

/*
                                      
          88                         88             
          ""                         88             
                                     88             
,adPPYba, 88 8b,dPPYba,   ,adPPYb,d8 88  ,adPPYba,  
I8[    "" 88 88P'   `"8a a8"    `Y88 88 a8P_____88  
 `"Y8ba,  88 88       88 8b       88 88 8PP"""""""  
aa    ]8I 88 88       88 "8a,   ,d88 88 "8b,   ,aa  
`"YbbdP"' 88 88       88  `"YbbdP"Y8 88  `"Ybbd8"'  
                          aa,    ,88                
                           "Y8bbdP"                 
                                    
                                    
                                    
                                    
8b,dPPYba,   ,adPPYba, 8b,dPPYba,   
88P'    "8a a8P_____88 88P'   `"8a  
88       d8 8PP""""""" 88       88  
88b,   ,a8" "8b,   ,aa 88       88  
88`YbbdP"'   `"Ybbd8"' 88       88  
88                                  
88

*/

//simulates the motion of a single pendulum
//using multiple finite difference methods
void single_pendulum()
{
	char processName[64] = "Single Pendulum";

	double t = 0;																//time
	double l = 9.81;														//length of pendulem in metres
	double m = 1.0;															//mass of pendulum in kg
	double gamma = 0.0;													//damping coefficient
	double beta = gamma/(m* sqrt( g*l ));				//matrix constant

	//initial conditions
	double initial_theta = 0.5;									//angle from vert (starting angle)
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

	std::cout.precision(1);											//sets the number of decimal places time is outputted to
																							//see: http://www.cplusplus.com/reference/ios/scientific/

	//set initial values
	euler_theta[0] = initial_theta;
	euler_w[0] = initial_w;

	leapfrog_theta[0] = initial_theta;
	leapfrog_w[0] = initial_w;
	
	rk4_theta[0] = initial_theta;
	rk4_w[0] = initial_w;



	//open up file and set column headings
	ofstream single_pen("data/single_pen_combined.csv");
	//first column is always time
	single_pen << "time" ;
	//euler
	single_pen << ",euler_theta,euler_w" ;
	//leapfrog
	single_pen << ",leapfrog_theta,leapfrog_w" ;
	//rk4
	single_pen << ",rk4_theta,rk4_w" ;
	//end column headers
	single_pen << "\n";

	//open up file and set column headings

	ofstream euler_single_pen("data/euler_single_pen.csv");
	ofstream leapfrog_single_pen("data/leapfrog_single_pen.csv");
	ofstream rk4_single_pen("data/rk4_single_pen.csv");

	//first column is always time
	euler_single_pen << "time" ;
	leapfrog_single_pen << "time" ;
	rk4_single_pen << "time" ;
	
	if(outputPositions)
	{
		//position column headings
		euler_single_pen << ",euler_theta,euler_w" ;
		leapfrog_single_pen << ",leapfrog_theta,leapfrog_w" ;
		rk4_single_pen << ",rk4_theta,rk4_w" ;
	}
	if(outputEnergies)
	{
		//energies column headings
		euler_single_pen << ",euler_T,euler_U, euler_E" ;
		leapfrog_single_pen << ",leapfrog_T,leapfrog_U, leapfrog_E" ;
		rk4_single_pen << ",rk4_T,rk4_U, rk4_E" ;
	}

	//end column headings
	euler_single_pen << "\n";
	leapfrog_single_pen << "\n";
	rk4_single_pen << "\n";

	for (int i = 0; i < numberOfSteps; i++)
	{
		/************** ENERGIES **************/

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

		/************** OUTPUT TO FILE **************/
		
		//output time
		t = h*i;
		euler_single_pen 		<< std::fixed << t ;
		leapfrog_single_pen << std::fixed << t ;
		rk4_single_pen 			<< std::fixed << t ;
		
		//positions
		if(outputPositions)
		{
			euler_single_pen << std::scientific << "," << euler_theta[i] << "," << euler_w[i] ;
			leapfrog_single_pen << std::scientific << "," << leapfrog_theta[i] << "," << leapfrog_w[i] ;
			rk4_single_pen << std::scientific << "," << rk4_theta[i] << "," << rk4_w[i] ;
		}
		//energies
		if(outputEnergies)
		{
			euler_single_pen << std::scientific << "," << euler_T[i] << "," << euler_U[i] << "," << euler_E[i];
			leapfrog_single_pen << std::scientific << "," << leapfrog_T[i] << "," << leapfrog_U[i] << "," << leapfrog_E[i];
			rk4_single_pen << std::scientific << "," << rk4_T[i] << "," << rk4_U[i] << "," << rk4_E[i];
		}
		euler_single_pen << "\n" ;
		leapfrog_single_pen << "\n" ;
		rk4_single_pen << "\n" ;

		//output time
		single_pen << std::fixed << t ;
		//output euler
		single_pen << std::scientific << "," << euler_theta[i] << "," << euler_w[i] ;
		//output leapfrog
		single_pen << std::scientific << "," << leapfrog_theta[i] << "," << leapfrog_w[i] ;
		//output rk4
		single_pen << std::scientific << "," << rk4_theta[i] << "," << rk4_w[i] ;
		//end output for this iteration
		single_pen << "\n";
		
		/************** EULER METHOD **************/

		//update euler_theta
		euler_theta[i+1] =  ( h*euler_w[i] ) + euler_theta[i];
		//update euler_w
		euler_w[i+1] = -h*euler_theta[i] + ( 1 - h*beta )*euler_w[i];

		/************** LEAPFROG METHOD **************/

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

		/************** RK4 METHOD **************/

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

		//update progress meter in terminal
		updateProgress(i, processName);

	}
} //end single_pendulum()

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

//simulates the motion of a double pendulum
//using the RK4 finite difference method
void double_pendulum()
{
	char processName[64] = "Double Pendulum";

	double t = 0;																//time (needed?)
	double l = 9.81;														//length of pendulem in metres
	double m = 1.0;															//mass of 1st pendulum in kg
	double M = 1.0;															//mass of 2nd pendulum in kg
	double gamma = 0.0;													//damping coefficient

	double beta = gamma/(m* sqrt( g*l ));				//matrix constant
	double G = beta;
	double R = M/m;															//matrix constant

	double initial_theta = 0.1;									//angle from vert for pendulum 1 (starting angle)
	double initial_psi = 0.0;										//angle from vert for pendulum 2 (starting angle)
	double initial_w = 0.0;											//rate of change of angle from vert for pendulum 1 (starting at rest)
	double initial_v = 0.0;											//rate of change of angle from vert for pendulum 2 (starting at rest)

	double rk4_theta [numberOfSteps];						//angular displacement of 1st pendulum
	double rk4_psi [numberOfSteps];							//angular displacement of 2nd pendulum
	double rk4_w [numberOfSteps];								//angular acceleration of 1st pendulum
	double rk4_v [numberOfSteps];								//angular acceleration of 2nd pendulum

	double y[4][numberOfSteps];									//hold all variables in multidimensional array

	//double k_1, k_2, k_3, k_4;								//used in old rk4 method
	double k1[4], k2[4], k3[4], k4[4];			//RK4 values

	std::cout.precision(1);											//sets the number of decimal places time is outputted to
																							//see: http://www.cplusplus.com/reference/ios/scientific/

	//set initial values
	rk4_theta[0] = initial_theta;
	rk4_psi[0] = initial_psi;
	rk4_w[0] = initial_w;
	rk4_v[0] = initial_v;

	//open up file and set column headings
	ofstream double_pen("data/double_pen.csv");
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

		/************** OLD RK4 METHOD ************** /

		//generate the k values for theta
		k_1 = h * ( rk4_w[i] );
		k_2 = h * ( rk4_w[i] + 0.5*k_1 );
		k_3 = h * ( rk4_w[i] + 0.5*k_2 );
		k_4 = h * ( rk4_w[i] + k_3);
		//generate next value of theta
		rk4_theta[i+1] = rk4_theta[i] + (1.0/6.0)*(k_1 + 2*k_2 + 2*k_3 + k_4);		//DONE

		//generate the k values for psi
		k_1 = h * ( rk4_v[i] );
		k_2 = h * ( rk4_v[i] + 0.5*k_1 );
		k_3 = h * ( rk4_v[i] + 0.5*k_2 );
		k_4 = h * ( rk4_v[i] + k_3);
		//generate next value of psi
		rk4_theta[i+1] = rk4_theta[i] + (1.0/6.0)*(k_1 + 2*k_2 + 2*k_3 + k_4);		//DONE

		//generate the k values for w
		k_1 = h * ( -1*(R+1)*rk4_theta[i] 					+ R*rk4_psi[i] 						- beta*rk4_w[i] );
		k_2 = h * ( -1*(R+1)*(rk4_theta[i] + 0.5*h) + R*(rk4_psi[i] + 0.5*h)  - beta*(rk4_w[i] + 0.5*h) );
		k_3 = h * ( -1*(R+1)*(rk4_theta[i] + 0.5*h) + R*(rk4_psi[i] + 0.5*h)  - beta*(rk4_w[i] + 0.5*h) );
		k_4 = h * ( -1*(R+1)*(rk4_theta[i] + h) 		+ R*(rk4_psi[i] + h) 			- beta*(rk4_w[i] + h) );
		//generate next value of w
		rk4_w[i+1] = rk4_w[i] + (1.0/6.0)*(k_1 + 2*k_2 + 2*k_3 + k_4);						//DONE
		
		//generate the k values for v
		k_1 = h * ( (R+1)*rk4_theta[i] 						- (R+1)*rk4_psi[i] 					 + beta*(1- (1/R))*rk4_w[i] 					- (beta/R)*rk4_v[i] );
		k_2 = h * ( (R+1)*(rk4_theta[i] + 0.5*h)  - (R+1)*(rk4_psi[i] + 0.5*h) + beta*(1- (1/R))*(rk4_w[i] + 0.5*h) - (beta/R)*(rk4_v[i] + 0.5*k_1) );
		k_3 = h * ( (R+1)*(rk4_theta[i] + 0.5*h)  - (R+1)*(rk4_psi[i] + 0.5*h) + beta*(1- (1/R))*(rk4_w[i] + 0.5*h) - (beta/R)*(rk4_v[i] + 0.5*k_2) );
		k_4 = h * ( (R+1)*(rk4_theta[i] + h) 			- (R+1)*(rk4_psi[i] + h) 		 + beta*(1- (1/R))*(rk4_w[i] + h) 		- (beta/R)*(rk4_v[i] + k_3) );
		//generate next value of v
		rk4_w[i+1] = rk4_w[i] + (1.0/6.0)*(k_1 + 2*k_2 + 2*k_3 + k_4);						//DONE

		/************** NEW RK4 METHOD **************/

		k1[0] = h * rk4_w[i];
		k1[1] = h * rk4_psi[i];
		k1[2] = h * ( -(R+1)*rk4_theta[i] + R*rk4_psi[i] - G*rk4_w[i] );
		k1[3] = h * ( (R+1)*rk4_theta[i] -(R+1)*rk4_psi[i] + G*(1-(1/R))*rk4_w[i] - (G/R)*rk4_v[i] );

		k2[0] = h * ( rk4_w[i] 		+ 0.5*k1[2] );
		k2[1] = h * ( rk4_psi[i]	+ 0.5*k1[3] );
		k2[2] = h * ( -(R+1)*(rk4_theta[i] + 0.5*k1[0]) + R*(rk4_psi[i] + 0.5*k1[1]) - G*(rk4_w[i] + 0.5*k1[2]) );
		k2[3] = h * ( (R+1)*(rk4_theta[i] + 0.5*k1[0]) -(R+1)*(rk4_psi[i] + 0.5*k1[1]) + G*(1-(1/R))*(rk4_w[i] + 0.5*k1[2]) - (G/R)*(rk4_v[i] + 0.5*k1[3]) );

		k3[0] = h * ( rk4_w[i] 		+ 0.5*k2[2] );
		k3[1] = h * ( rk4_psi[i]	+ 0.5*k2[3] );
		k3[2] = h * ( -(R+1)*(rk4_theta[i] + 0.5*k2[0]) + R*(rk4_psi[i] + 0.5*k2[1]) - G*(rk4_w[i] + 0.5*k2[2]) );
		k3[3] = h * ( (R+1)*(rk4_theta[i] + 0.5*k2[0]) -(R+1)*(rk4_psi[i] + 0.5*k2[1]) + G*(1-(1/R))*(rk4_w[i] + 0.5*k2[2]) - (G/R)*(rk4_v[i] + 0.5*k1[3]) );

		k4[0] = h * ( rk4_w[i] 		+ k3[2] );
		k4[1] = h * ( rk4_psi[i]	+ k3[3] );
		k4[2] = h * ( -(R+1)*(rk4_theta[i] + k3[0]) + R*(rk4_psi[i] + k3[1]) - G*(rk4_w[i] + k3[2]) );
		k4[3] = h * ( (R+1)*(rk4_theta[i] + k3[0]) -(R+1)*(rk4_psi[i] + k3[1]) + G*(1-(1/R))*(rk4_w[i] + k3[2]) - (G/R)*(rk4_v[i] + k3[3]) );

		rk4_theta[i+1]	= rk4_theta[i]	+ (1.0/6.0)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
		rk4_psi[i+1]		= rk4_psi[i]		+ (1.0/6.0)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
		rk4_w[i+1]			= rk4_w[i]			+ (1.0/6.0)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
		rk4_v[i+1]			= rk4_v[i]			+ (1.0/6.0)*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
		
		/************** PROGRESS **************/
		updateProgress(i, processName);

	}
} //end double_pendulum

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