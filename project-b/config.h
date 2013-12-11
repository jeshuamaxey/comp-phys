/*********************************************
*	SELECT SIMULATIONS TO RUN
**********************************************/

#define simulateForZeroB true
#define simulateForNonZeroB false
#define simulateForVaryingB false

/*********************************************
*	PROGRAM PARAMETERS
**********************************************/

//Size of single dimension of mesh
#define N 20
//Number of previous energies kept track of to test for equilibrium
#define numberPrevEs 3*N
//number seeds used from rng
#define numberOfSeeds 10
//The threshold for percent difference between moving averages.
//If a difference is below this value, the system is said to have reached equilibrium
#define pc_diff_max 0.0001
//number of times N*N flips are posited after equilibrium is found
#define simulationsPastEquilibrium N*N

//J/(k_b*T) = beta
#define J 1.0

/*********************************************
*	FOR VARYING BETA
**********************************************/

//Increment in beta between simulations
#define beta_step 0.005
//Minimum value of beta simulated
#define beta_min 0.15
//Maximum value of beta simulated
#define beta_max 1.0

/*********************************************
*	FOR VARYING mu_B
**********************************************/

//Increment in beta between simulations
#define mu_B_step 1.0*J
//Minimum value of beta simulated
#define mu_B_min -0.5*J
//Maximum value of beta simulated
#define mu_B_max 0.5*J

/*********************************************
*	DETERMINE WHAT THE PROGRAM OUTPUTS
**********************************************/
#define outputE false
#define outputM true
#define outputSHC false
#define outputMS false

#define outputOneOverBeta true //used to affect the type of plot