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

extern "C" struct OpesciProfiling {
${define_profiling}
};

extern "C" int opesci_execute(OpesciGrid *grid, OpesciProfiling *profiling) {
% if profiling==True:
int err = opesci_papi_init();
% endif

${define_constants}
${declare_fields}

#pragma omp parallel
{
% if profiling==True:
float real_time;
float proc_time;
float mflops;
long long flpins;
opesci_flops(&real_time, &proc_time, &flpins, &mflops);
% endif

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

% if profiling==True:
opesci_flops(&real_time, &proc_time, &flpins, &mflops);
#pragma omp critical
{
profiling->g_rtime = fmax(profiling->g_rtime, real_time);
profiling->g_ptime = fmax(profiling->g_ptime, proc_time);
profiling->g_mflops += mflops;
}
% endif

} // end of parallel section

${store_fields}

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
OpesciProfiling profiling;

 opesci_execute(&grid, &profiling);
opesci_convergence(&grid, &conv);
${print_convergence}

% if profiling==True:
printf("PAPI:: Max real_time: %f (sec)\n", profiling.g_rtime);
printf("PAPI:: Max proc_time: %f (sec)\n", profiling.g_ptime);
printf("PAPI:: Total MFlops/s: %f\n", profiling.g_mflops);
% endif

return 0;
}
