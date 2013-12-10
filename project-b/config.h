/*********************************************
*	PROGRAM PARAMETERS
**********************************************/

//Size of single dimension of mesh
#define N 100
//Increment in beta between simulations
#define beta_step 0.01
//Minimum value of beta simulated
#define beta_min 1.0
//Maximum value of beta simulated
#define beta_max 1.0
//Number of previous energies kept track of to test for equilibrium
#define numberPrevEs 500
//The threshold for percent difference between moving averages.
//If a difference is below this value, the system is said to have reached equilibrium
#define pc_diff_max 0.0001

/*********************************************
* DETERMINE WHAT THE PROGRAM OUTPUTS
**********************************************/
#define outputE false
#define outputM false
#define outputSHC false
#define outputMS true