<%include file="copyright.txt"/>

<%include file="common_include.txt"/>
#include <cmath>

int main(){

  // defined constants
  float rho = 1.0;
  float mu = 0.5;
  float lambda = 0.5;
  float _tmax = 1.0;
  float h = 0.025;

  // calculated constatns
  float Vp = sqrt((lambda+2*mu)/rho);
  float Vs = sqrt(mu/rho);
  float dt = 0.5*h/Vp; // need to check
  float beta = 1.0/rho;
  int ntsteps = (int) _tmax/dt;
  int margin = 2; // ghost cells for boundary conditions
  int dimx, dimy;
  dimx = dimy = (int) 1.0/h + 1 + margin*2;

  // time periodicity for update
  const int _tp = ${time_period};

  // set up solution mesh
  float U[_tp][dimx][dimy];
  float V[_tp][dimx][dimy];
  float Txx[_tp][dimx][dimy];
  float Tyy[_tp][dimx][dimy];
  float Txy[_tp][dimx][dimy];

  // shared variables
  int t, t1;

#pragma omp parallel
  {
  // Initialise fields
  ${initialise}

  } // end of parallel section

  return 0;
}
