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
% if pluto==True:
<%include file ="pluto_include.txt"/>
% endif

extern "C" struct OpesciGrid {
${define_fields}
};

extern "C" struct OpesciConvergence {

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
% if numevents_papi>0:
${define_papi_events}
opesci_papi_start_counters(numevents, events);
% else:
float real_time;
float proc_time;
float mflops;
long long flpins;
opesci_flops(&real_time, &proc_time, &flpins, &mflops);
% endif
% endif

${initialise}

for(int _ti=0;_ti<ntsteps;_ti++){

${time_stepping}

% if pluto==True:
{
#pragma scop
% endif
${primary_loop}
% if pluto==True:
#pragma endscop
}
% endif

} // end of time loop

% if profiling==True:
% if numevents_papi>0:
opesci_papi_read_counters(numevents, counters);
#pragma omp critical
{
${sum_papi_events}
}
% else:
opesci_flops(&real_time, &proc_time, &flpins, &mflops);
#pragma omp critical
{
profiling->g_rtime = fmax(profiling->g_rtime, real_time);
profiling->g_ptime = fmax(profiling->g_ptime, proc_time);
profiling->g_mflops += mflops;
}
% endif
% endif

} // end of parallel section

${store_fields}

return 0;
}

extern "C" int opesci_convergence(OpesciGrid *grid, OpesciConvergence *conv) {
${define_constants}
${load_fields}

return 0;
}

int main(){
OpesciGrid grid;
OpesciConvergence conv;
OpesciProfiling profiling;

 opesci_execute(&grid, &profiling);


% if profiling==True:
printf("PAPI:: Max real_time: %f (sec)\n", profiling.g_rtime);
printf("PAPI:: Max proc_time: %f (sec)\n", profiling.g_ptime);
printf("PAPI:: Total MFlops/s: %f\n", profiling.g_mflops);
% endif

return 0;
}
