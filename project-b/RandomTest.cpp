#include <iostream>
#include <iomanip>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

int main (void) {
  const gsl_rng_type * T;
  gsl_rng * r;

  int i, n = 10;

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);


  std::cout << "Default seed" << std::endl;
  for (i = 0; i < n; i++) {
  	double u = gsl_rng_uniform (r);
  	std::cout << std::setw(10) << u << std::endl;
  }

  gsl_rng_set(r,0);
  std::cout << "seed=0 (Should repeat same sequence)" << std::endl;
  for (i = 0; i < n; i++) {
  	double u = gsl_rng_uniform (r);
  	std::cout << std::setw(10) << u << std::endl;
  }

  gsl_rng_set(r,(unsigned long) 116426264);
  std::cout << "seed = 11642626462 (Should create new sequence)" 
			<< std::endl;
  for (i = 0; i < n; i++) {
  	double u = gsl_rng_uniform (r);
  	std::cout << std::setw(10) << u << std::endl;
  }

  return 0;
}
