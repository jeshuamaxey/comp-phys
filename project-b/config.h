/*********************************************
*	SELECT SIMULATIONS TO RUN
**********************************************/

#define simulateForZeroB false							// all of these plots are as expected :)
#define simulateForNonZeroB false
#define simulateForVaryingB true

/*********************************************
*	PROGRAM PARAMETERS
**********************************************/

//Size of single dimension of mesh
#define N 50
//Number of previous energies kept track of to test for equilibrium
#define numberPrevEs 3*N
//number seeds used from rng (max: 15)
#define numberOfSeeds 15
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
#define mu_B_step 0.005
//Minimum value of beta simulated
#define mu_B_min -1
//Maximum value of beta simulated
#define mu_B_max 1
//set the fixed value of beta for this simulation
#define fixedBeta 0.41

/*********************************************
*	DETERMINE WHAT THE PROGRAM OUTPUTS
**********************************************/
#define outputE true
#define outputM false
#define outputSHC false
#define outputMS false

#define outputOneOverBeta false //used to affect the type of plot


/*
while (i <= limit) {
		//set magnetic field strength
		if(mu_B_inRangeOfInterest(mu_B) && i !=0.0 )
		{
			mu_B = mu_B_min+(i*zoom*mu_B_step);
			i+=zoom;
		} else
		{
			mu_B = mu_B_min+(i*mu_B_step);
			i++;
		}
		//run simulation under these conditions
		runSimulation(beta, mu_B, outputFile1, a, 2);
		//
		updateProgress(i/(mu_B_max*2/mu_B_step), a);
	}



*/
