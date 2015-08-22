<%include file="copyright.txt"/>
#ifdef _MSC_VER
#define M_PI 3.14159265358979323846
#endif
<%include file="common_include.txt"/>
% if io==True:
<%include file="io_include.txt"/>
% endif
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