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
//void outputToFile(ostream&);								//not used
void single_pendulum();
void double_pendulum();
void updateProgress(int);
void updateProgress(int, char[]);

int main() {
	single_pendulum();
	double_pendulum();
	/************** PROGRESS METER **************/
	std::cout << "\r"<< "100\% - all done. Simulated time elapsed: " << simulatedTime << "s\n" << std::flush;
	return 0;
}

void single_pendulum()
{
	char processName[64] = "Single Pendulum";

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

		updateProgress(i, processName);

	}
} //end single_pendulum()

void double_pendulum()
{
	char processName[64] = "Double Pendulum";

	double t = 0;																//time (needed?)
	double l = 9.81;														//length of pendulem in metres
	double m = 1.0;															//mass of 1st pendulum in kg
	double M = 1.0;															//mass of 2nd pendulum in kg
	double gamma = 1.0;													//damping coefficient

	double beta = gamma/(m* sqrt( g*l ));				//matrix constant
	double R = M/m;															//matrix constant

	double initial_theta = 0.2;									//angle from vert for pendulum 1 (starting angle)
	double initial_psi = 0.2;										//angle from vert for pendulum 2 (starting angle)
	double initial_w = 0.0;											//rate of change of angle from vert for pendulum 1 (starting at rest)
	double initial_v = 0.0;											//rate of change of angle from vert for pendulum 2 (starting at rest)

	double rk4_theta [numberOfSteps];						//angular displacement of 1st pendulum
	double rk4_psi [numberOfSteps];							//angular displacement of 2nd pendulum
	double rk4_w [numberOfSteps];								//angular acceleration of 1st pendulum
	double rk4_v [numberOfSteps];								//angular acceleration of 2nd pendulum

	double k_1, k_2, k_3, k_4;									//RK4 values

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

		/************** RK4 METHOD **************/

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
		rk4_w[i+1] = rk4_w[i] + (1.0/6.0)*(k_1 + 2*k_2 + 2*k_3 + k_4);
		
		//generate the k values for v
		k_1 = h * ( (R+1)*rk4_theta[i] 						- (R+1)*rk4_psi[i] 					 + beta*(1- (1/R))*rk4_w[i] 					- (beta/R)*rk4_v[i] );
		k_2 = h * ( (R+1)*(rk4_theta[i] + 0.5*h)  - (R+1)*(rk4_psi[i] + 0.5*h) + beta*(1- (1/R))*(rk4_w[i] + 0.5*h) - (beta/R)*(rk4_v[i] + 0.5*k_1) );
		k_3 = h * ( (R+1)*(rk4_theta[i] + 0.5*h)  - (R+1)*(rk4_psi[i] + 0.5*h) + beta*(1- (1/R))*(rk4_w[i] + 0.5*h) - (beta/R)*(rk4_v[i] + 0.5*k_2) );
		k_4 = h * ( (R+1)*(rk4_theta[i] + h) 			- (R+1)*(rk4_psi[i] + h) 		 + beta*(1- (1/R))*(rk4_w[i] + h) 		- (beta/R)*(rk4_v[i] + k_3) );
		//generate next value of v
		rk4_w[i+1] = rk4_w[i] + (1.0/6.0)*(k_1 + 2*k_2 + 2*k_3 + k_4);						//DONE
		
		/************** PROGRESS **************/
		updateProgress(i, processName);

	}
} //end double_pendulum

//updates the progress display on commandline
void updateProgress(int i, char *name) {
		cout << std::fixed << "\r"<< "Currently running " << name << " - " << (float(i)/numberOfSteps)*100 << "\%" << flush;
}