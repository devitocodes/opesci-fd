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
  #pragma omp for
  for(int x=margin;x<dimx-margin;x++){
    for(int y=margin;y<dimy-margin;y++){
      float xx = (x-margin)*h;
      float yy = (y-margin)*h;
      ${Txx_init};
      ${Tyy_init};
    }
  }
  #pragma omp for
  for(int x=margin;x<dimx-margin-1;x++){
    for(int y=margin;y<dimy-margin;y++){
      float xx = (x-margin+0.5)*h;
      float yy = (y-margin)*h;
      ${U_init};
    }
  }
  #pragma omp for
  for(int x=margin;x<dimx-margin;x++){
    for(int y=margin;y<dimy-margin-1;y++){
      float xx = (x-margin)*h;
      float yy = (y-margin+0.5)*h;
      ${V_init};
    }
  }
  #pragma omp for
  for(int x=margin;x<dimx-margin-1;x++){
    for(int y=margin;y<dimy-margin-1;y++){
      float xx = (x-margin+0.5)*h;
      float yy = (y-margin+0.5)*h;
      ${Txy_init};
    }
  }

    // main time loop
  for(int _ti=0;_ti<ntsteps;_ti++){

    // shared variables
    #pragma omp single
    {
      t = _ti % _tp; // array index of current time step
      t1 = (t+1) % _tp; // array index of the grid to be updated
    }

    // Compute stresses
    #pragma omp for
    for(int x=margin;x<dimx-margin;x++){
      for(int y=margin;y<dimy-margin;y++){
        ${Txx};
        ${Tyy};
        ${Txy};
      }
    }

    // update ghost cells for boundary conditions
    #pragma omp for
    for(int x=0;x<dimx;x++){
      // boundary y=2
      Tyy[t1][x][1] = -Tyy[t1][x][3];
      Txy[t1][x][0] = -Txy[t1][x][3];
      Txy[t1][x][1] = -Txy[t1][x][2];
      // boundary y=dimy+2
      Tyy[t1][x][dimy-2] = -Tyy[t1][x][dimy-4];
      Txy[t1][x][dimy-1] = -Txy[t1][x][dimy-4];
      Txy[t1][x][dimy-2] = -Txy[t1][x][dimy-3];
    }
    #pragma omp for
    for(int y=0;y<dimy;y++){
      // boundary x=2
      Tyy[t1][1][y] = -Tyy[t1][3][y];
      Txy[t1][0][y] = -Txy[t1][3][y];
      Txy[t1][1][y] = -Txy[t1][2][y];
      // boundary x=dimx+2
      Tyy[t1][dimx-2][y] = -Tyy[t1][dimx-4][y];
      Txy[t1][dimx-1][y] = -Txy[t1][dimx-4][y];
      Txy[t1][dimx-2][y] = -Txy[t1][dimx-3][y];
    }

    // Compute velocities
    #pragma omp for
    for(int x=margin;x<dimx-margin;x++){
      for(int y=margin;y<dimy-margin;y++){
        ${U};
        ${V};
      }
    }

    // update ghost cells for boundary conditions
    #pragma omp for
    for(int x=0;x<dimx;x++){
      // boundary y=2
      ${bc_U_y0};
      ${bc_V_y0};
      // boundary y=dimy+2
      ${bc_U_y1};
      ${bc_V_y1};
    }
    #pragma omp for
    for(int y=0;y<dimy;y++){
      // boundary x=2
      ${bc_U_x0};
      ${bc_V_x0};
      // boundary x=dimy-3
      ${bc_U_x1};
      ${bc_V_x1};
    }
  } // end of time loop
  } // end of parallel section

  return 0;
}
