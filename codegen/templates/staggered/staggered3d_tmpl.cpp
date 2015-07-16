<%include file="copyright.txt"/>

<%include file="common_include.txt"/>
#include <cmath>
#include <cstdio>
#include <string>

int main(){

  ${define_constants}
  ${declare_fields}

#pragma omp parallel
  {
  ${initialise}
  ${initialise_bc}

  for(int _ti=0;_ti<ntsteps;_ti++){
    
    ${time_stepping}

    ${stress_loop}
    ${stress_bc}
    ${velocity_loop}
    ${velocity_bc}

    ${output_step}

  } // end of time loop
  } // end of parallel section

  ${output_final}

  return 0;
}