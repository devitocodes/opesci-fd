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

  float Txx_diff = 0.0;
  float Tyy_diff = 0.0;
  float Txy_diff = 0.0;
  float U_diff = 0.0;
  float V_diff = 0.0;
  float tf1 = ntsteps*dt;
  float tf2 = tf1 + dt/2;

  ${converge_test}

  std::cout<<Txx_diff<<'\n';
  std::cout<<Tyy_diff<<'\n';
  std::cout<<Txy_diff<<'\n';
  std::cout<<U_diff<<'\n';
  std::cout<<V_diff<<'\n';
  return 0;
}
