#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#define simulatedTime 12.0											//simulated time (seconds)
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

//helper functions
void updateProgress(int);
void updateProgress(int, char[]);
void done();

int main()
{
	single_pendulum();
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