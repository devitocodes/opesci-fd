<%include file="copyright.txt"/>
#ifdef _MSC_VER
#define M_PI 3.14159265358979323846
#endif
<%include file="common_include.txt"/>
% if io==True:
<%include file="io_include.txt"/>
% endif
% if profiling==True:
#include "opesciProfiling.h"
float real_time;
float proc_time;
float mflops;
long long flpins;
% endif

#include <cmath>
#include <cstdio>
#include <string>

extern "C" struct OpesciGrid {
${define_fields}
};

extern "C" struct OpesciConvergence {
${define_convergence}
};

extern "C" int opesci_execute(OpesciGrid *grid) {

${define_constants}
${declare_fields}

% if profiling==True:
opesci_flops(&real_time, &proc_time, &flpins, &mflops);
% endif

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

${store_fields}

% if profiling==True:
opesci_flops(&real_time, &proc_time, &flpins, &mflops);
printf("PAPI:: Total Flops:\n");
printf("PAPI:: Total rtime: %f (sec)\n", (real_time));
printf("PAPI:: Total ptime: %f (sec)\n", (proc_time));
printf("PAPI:: MFlops/s: %f\n", mflops);
% endif

return 0;
}

extern "C" int opesci_convergence(OpesciGrid *grid, OpesciConvergence *conv) {
${define_constants}
${load_fields}

${converge_test}
return 0;
}

int main(){
OpesciGrid grid;
OpesciConvergence conv;

opesci_execute(&grid);
opesci_convergence(&grid, &conv);
${print_convergence}

return 0;
}
