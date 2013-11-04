#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#define simulatedTime 30											//simulated time (seconds)
#define h 0.005																//step size
#define numberOfSteps  int(simulatedTime / h)	//used for sizing arrays
#define g 9.81																//acceleration due to gravity
#define pi atan(1.0)													//mutha fuckin pi man

using namespace std;

//function prototypes

int main()
{
	double t = 0;																//time (needed?)
	double l = 9.81;														//length of pendulem in metres
	double m = 1.0;															//mass of pendulum in kg
	double gamma = 1.0;													//damping coefficient
	double beta = gamma/(m* sqrt( g*l ));				//matrix constant

	double initial_theta = 0.5;									//angle from vert (starting angle)
	double initial_w = 0;												//rate of change of angle from vert (starting at rest)
	
	double euler_theta_n, euler_w_n;						//used to keep track of current values during update calculations
	double euler_theta [numberOfSteps];					//stores all theta values
	double euler_w [numberOfSteps];							//stores all w values

	double leapfrog_theta_n, leapfrog_w_n;			//used to keep track of current values during update calculations
	double leapfrog_theta [numberOfSteps];			//stores all theta values
	double leapfrog_w [numberOfSteps];					//stores all w values

	double rk4_theta [numberOfSteps];
	double rk4_w [numberOfSteps];
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
	ofstream single_pen("data/single_pen.csv");
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

	for (int i = 0; i < numberOfSteps; i++)
	{
		/************** OUTPUT **************/
		//output time
		single_pen << std::fixed << h*i ;
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

		/************** PROGRESS METER **************/
		cout << std::fixed<< "\r"<< (float(i)/numberOfSteps)*100 << "\%" << flush;

	}

	/************** PROGRESS METER **************/
	std::cout << "\r"<< "100\% - all done. Simulated time elapsed: " << numberOfSteps*h << "s\n" << std::flush;
	return 0;
}