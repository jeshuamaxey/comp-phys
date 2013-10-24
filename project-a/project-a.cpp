#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#define numberOfSteps 5000									//number of time steps
#define h 0.005														//step size
#define g 9.81														//acceleration due to gravity

using namespace std;

//function prototypes

int main()
{
	double t = 0;																//time (needed?)
	double l = 1.0;															//length of pendulem in metres
	double m = 1.0;															//mass of pendulum in kg
	double gamma = 200.0;												//damping coefficient
	double beta = gamma/(m* sqrt( g*l ));				//matrix constant

	double initial_theta = 0.3;									//angle from vert (starting angle)
	double initial_w = 0;												//rate of change of angle from vert (starting at rest)
	
	double euler_theta_n, euler_w_n;						//used to keep track of current values during update calculations
	double euler_theta [numberOfSteps];					//stores all theta values
	double euler_w [numberOfSteps];							//stores all w values

	double leapfrog_theta_n, leapfrog_w_n;			//used to keep track of current values during update calculations
	double leapfrog_theta [numberOfSteps];			//stores all theta values
	double leapfrog_w [numberOfSteps];					//stores all w values

	
	//set initial values
	euler_theta[0] = initial_theta;
	euler_w[0] = initial_w;
	leapfrog_theta[0] = initial_theta;
	leapfrog_w[0] = initial_w;

	//open up file and set column headings
	ofstream single_pen("data/single_pen.csv");
	single_pen << "time,euler_theta,leapfrog_theta" << endl; //euler only for now
	//single_pen << "time,euler_theta,euler_w,leapfrog_theta,leapfrog_w" << endl;

	for (int i = 0; i < numberOfSteps; i++)
	{
		/************** OUTPUT **************/

		//output results to file
		//single_pen << i*h << "," << euler_theta[i] << "," << euler_w[i] << "," << leapfrog_theta[i] << "," << leapfrog_w[i] << endl;
		single_pen << i*h << "," << euler_theta[i] << "," << leapfrog_theta[i] << endl;
		/************** EULER METHOD **************/

		//update euler_theta
		euler_theta[i+1] =  ( h*euler_w[i] ) + euler_theta[i];
		//update euler_w
		euler_w[i+1] = ( (1-h) * euler_w[i] ) - ( h * beta * euler_theta[i] );

		/************** LEAPFROG METHOD ************** /

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
			leapfrog_theta[i+1] =  leapfrog_theta[i-1] + leapfrog_w[i]*2*h;
			//update leapfrog_w
			leapfrog_w[i+1] = leapfrog_w[i-1] - ( leapfrog_w[i] + ( beta*leapfrog_theta[i] )*2*h );
		}

		if(i>400)
		{
			break;
		}

		/************** PROGRESS METER **************/
		cout << "\r"<< (float(i)/numberOfSteps)*100 << "\%" << flush;

	}

	/************** PROGRESS METER **************/
	std::cout << "\r"<< "100\% - all done. Simulated time elapsed: " << numberOfSteps*h << "s\n" << std::flush;
	return 0;
}