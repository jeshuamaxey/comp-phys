/*********************************************
*	PROGRAM PARAMETERS
**********************************************/

//Size of single dimension of mesh
#define N 10
//Increment in beta between simulations
#define beta_step 0.005
//Minimum value of beta simulated
#define beta_min 0.15
//Maximum value of beta simulated
#define beta_max 1.0
//Number of previous energies kept track of to test for equilibrium
#define numberPrevEs 3*N
//The threshold for percent difference between moving averages.
//If a difference is below this value, the system is said to have reached equilibrium
#define pc_diff_max 0.0001
//
#define simulationsPastEquilibrium 500

/*********************************************
*	DETERMINE WHAT THE PROGRAM OUTPUTS
**********************************************/
#define outputE false
#define outputM false
#define outputSHC true	//CHECK it
#define outputMS false

#define outputOneOverBeta true //used to affect the type of plot