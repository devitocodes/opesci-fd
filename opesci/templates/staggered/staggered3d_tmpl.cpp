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

extern "C" struct OpesciGrid {
${define_fields}
};

extern "C" int opesci_execute(OpesciGrid *grid) {

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

${store_fields}

return 0;
}

extern "C" int opesci_convergence(OpesciGrid *grid) {
${define_constants}
${load_fields}

${converge_test}
}

int main(){
  OpesciGrid grid;

  opesci_execute(&grid);
  opesci_convergence(&grid);

  return 0;
}
