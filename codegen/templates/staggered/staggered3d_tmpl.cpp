<%include file="copyright.txt"/>

<%include file="common_include.txt"/>
#include <cmath>
#include <cstdio>

int main(){

  const int _tp = 2;
  // defined constants
  ${define_constants}
  // set up solution mesh
  ${declare_fields}

#pragma omp parallel
  {
  // Initialise fields
  ${initialise}
  ${initialise_bc}

  for(int _ti=0;_ti<ntsteps;_ti++){
    
    // shared variables
    #pragma omp single
    {
      t = _ti % _tp; // array index of current time step
      t1 = (t+1) % _tp; // array index of the grid to be updated
    }

    ${stress_loop}

    ${stress_bc}

    ${velocity_loop}

    ${velocity_bc}

  } // end of time loop
  } // end of parallel section

  ${output}

  return 0;
}
