/*********************************************
*	SELECT SIMULATIONS TO RUN
**********************************************/

#define simulateForZeroB false							// all of these plots are as expected :)
#define simulateForNonZeroB true
#define simulateForVaryingB false

/*********************************************
*	PROGRAM PARAMETERS
**********************************************/

//Size of single dimension of mesh
#define N 20
//Number of previous energies kept track of to test for equilibrium
#define numberPrevEs 3*N
//number seeds used from rng (max: 15)
#define numberOfSeeds 10
//The threshold for percent difference between moving averages.
//If a difference is below this value, the system is said to have reached equilibrium
#define pc_diff_max 0.0001
//number of times N*N flips are posited after equilibrium is found
#define simulationsPastEquilibrium N*N

//J/(k_b*T) = beta
#define J 1.0 				//not used in code any longer

/*********************************************
*	FOR VARYING BETA
**********************************************/

//Increment in beta between simulations
#define beta_step 0.005
//Minimum value of beta simulated
#define beta_min 0.001
//Maximum value of beta simulated
#define beta_max 0.6
//
#define nonZero_mu_B 0.5

/*********************************************
*	FOR VARYING mu_B
**********************************************/

//Increment in beta between simulations
#define mu_B_step 0.1
//Minimum value of beta simulated
#define mu_B_min -1
//Maximum value of beta simulated
#define mu_B_max 1
//set the fixed value of beta for this simulation
#define fixedBeta 0.01

/*********************************************
*	DETERMINE WHAT THE PROGRAM OUTPUTS
**********************************************/
#define outputE false
#define outputM false
#define outputSHC false
#define outputMS true

#define outputOneOverBeta false //used to affect the type of plot
