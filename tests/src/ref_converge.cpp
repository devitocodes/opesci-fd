/*  Copyright (C) 2015 Imperial College London and others.
 *
 *  Please see the AUTHORS file in the main source directory for a
 *  full list of copyright holders.
 *
 *  Gerard Gorman
 *  Department of Earth Science and Engineering
 *  Imperial College London
 *
 *  g.gorman@imperial.ac.uk
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer in the documentation and/or other materials provided
 *  with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 *  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 *  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 *  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 *  SUCH DAMAGE.
 */

/* Test case implements:

   Graves, Robert W. "Simulating seismic wave propagation in 3D
   elastic media using staggered-grid finite differences." Bulletin of
   the Seismological Society of America 86.4 (1996): 1091-1106.

*/


//#include "opesciIO.h"
//#include "opesciHandy.h"

#include <cassert>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
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
  const int _tp = 2;

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
      Txx[0][x][y]=0;
      Tyy[0][x][y]=0;
    }
  }
  #pragma omp for
  for(int x=margin;x<dimx-margin-1;x++){
    for(int y=margin;y<dimy-margin;y++){
      float xx = (x-margin+0.5)*h;
      float yy = (y-margin)*h;
      U[0][x][y]=1.4142135623731*M_PI*Vs*sin(M_PI*yy)*cos(M_PI*xx)*cos(0.707106781186548*M_PI*Vs*dt);
    }
  }
  #pragma omp for
  for(int x=margin;x<dimx-margin;x++){
    for(int y=margin;y<dimy-margin-1;y++){
      float xx = (x-margin)*h;
      float yy = (y-margin+0.5)*h;
      V[0][x][y]=-1.4142135623731*M_PI*Vs*sin(M_PI*xx)*cos(M_PI*yy)*cos(0.707106781186548*M_PI*Vs*dt);
    }
  }
  #pragma omp for
  for(int x=margin;x<dimx-margin-1;x++){
    for(int y=margin;y<dimy-margin-1;y++){
      float xx = (x-margin+0.5)*h;
      float yy = (y-margin+0.5)*h;
      Txy[0][x][y]=0.0;
    }
  }

  // initialise boundary
  #pragma omp for
  for(int x=0;x<dimx;x++){
    // boundary y=2
    Tyy[0][x][1] = -Tyy[0][x][3];
    Txy[0][x][0] = -Txy[0][x][3];
    Txy[0][x][1] = -Txy[0][x][2];
    // boundary y=dimy+2
    Tyy[0][x][dimy-2] = -Tyy[0][x][dimy-4];
    Txy[0][x][dimy-1] = -Txy[0][x][dimy-4];
    Txy[0][x][dimy-2] = -Txy[0][x][dimy-3];
  }
  #pragma omp for
  for(int y=0;y<dimy;y++){
    // boundary x=2
    Tyy[0][1][y] = -Tyy[0][3][y];
    Txy[0][0][y] = -Txy[0][3][y];
    Txy[0][1][y] = -Txy[0][2][y];
    // boundary x=dimx+2
    Tyy[0][dimx-2][y] = -Tyy[0][dimx-4][y];
    Txy[0][dimx-1][y] = -Txy[0][dimx-4][y];
    Txy[0][dimx-2][y] = -Txy[0][dimx-3][y];
  }
  #pragma omp for
  for(int x=0;x<dimx;x++){
    // boundary y=2
    U[0][x][1]=U[0][x][3] - V[0][x][1] - V[0][x][2] + V[0][x + 1][1] + V[0][x + 1][2];
    V[0][x][1]=(lambda*U[0][x][2] - lambda*U[0][x - 1][2] + lambda*V[0][x][2] + 2*mu*V[0][x][2])/(lambda + 2*mu);
    // boundary y=dimy+2
    U[0][x][dimy - 2]=U[0][x][dimy - 4] + V[0][x][dimy - 4] + V[0][x][dimy - 3] - V[0][x + 1][dimy - 4] - V[0][x + 1][dimy - 3];
    V[0][x][dimy - 3]=(-lambda*U[0][x][dimy - 3] + lambda*U[0][x - 1][dimy - 3] + lambda*V[0][x][dimy - 4] + 2*mu*V[0][x][dimy - 4])/(lambda + 2*mu);
  }
  #pragma omp for
  for(int y=0;y<dimy;y++){
    // boundary x=2
    U[0][1][y]=(lambda*(U[0][2][y] + V[0][2][y] - V[0][2][y - 1]) + 2*mu*V[0][2][y] - 2*mu*V[0][2][y - 1])/lambda;
    V[0][1][y]=-U[0][1][y] + U[0][1][y + 1] - U[0][2][y] + U[0][2][y + 1] + V[0][3][y];
    // boundary x=dimy-3
    U[0][dimx - 3][y]=(lambda*(U[0][dimx - 4][y] - V[0][dimx - 3][y] + V[0][dimx - 3][y - 1]) - 2*mu*V[0][dimx - 3][y] + 2*mu*V[0][dimx - 3][y - 1])/lambda;
    V[0][dimx - 2][y]=U[0][dimx - 4][y] - U[0][dimx - 4][y + 1] + U[0][dimx - 3][y] - U[0][dimx - 3][y + 1] + V[0][dimx - 4][y];
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
    for(int x=margin+1;x<dimx-margin-1;x++){
      for(int y=margin+1;y<dimy-margin-1;y++){
        Txx[t1][x][y]=(1.0F/24.0F)*dt*(27*lambda*U[t][x][y] + lambda*U[t][x - 2][y] - 27*lambda*U[t][x - 1][y] - lambda*U[t][x + 1][y] + 27*lambda*V[t][x][y] + lambda*V[t][x][y - 2] - 27*lambda*V[t][x][y - 1] - lambda*V[t][x][y + 1] + 54*mu*U[t][x][y] + 2*mu*U[t][x - 2][y] - 54*mu*U[t][x - 1][y] - 2*mu*U[t][x + 1][y])/h + Txx[t][x][y];
        Tyy[t1][x][y]=(1.0F/24.0F)*dt*(27*lambda*U[t][x][y] + lambda*U[t][x - 2][y] - 27*lambda*U[t][x - 1][y] - lambda*U[t][x + 1][y] + 27*lambda*V[t][x][y] + lambda*V[t][x][y - 2] - 27*lambda*V[t][x][y - 1] - lambda*V[t][x][y + 1] + 54*mu*V[t][x][y] + 2*mu*V[t][x][y - 2] - 54*mu*V[t][x][y - 1] - 2*mu*V[t][x][y + 1])/h + Tyy[t][x][y];
        Txy[t1][x][y]=(1.0F/24.0F)*dt*mu*(-27*U[t][x][y] + U[t][x][y - 1] + 27*U[t][x][y + 1] - U[t][x][y + 2] - 27*V[t][x][y] + V[t][x - 1][y] + 27*V[t][x + 1][y] - V[t][x + 2][y])/h + Txy[t][x][y];
      }
    }

    // update ghost cells for boundary conditions
    #pragma omp for
    for(int x=0;x<dimx;x++){
      // boundary y=2
      Txx[t1][x][1] = -Txx[t1][x][3];
      Txx[t1][x][2] = 0.0;
      Tyy[t1][x][1] = -Tyy[t1][x][3];
      Tyy[t1][x][2] = 0.0;
      Txy[t1][x][0] = -Txy[t1][x][3];
      Txy[t1][x][1] = -Txy[t1][x][2];
      // boundary y=dimy+2
      Txx[t1][x][dimy-2] = -Txx[t1][x][dimy-4];
      Txx[t1][x][dimy-3] = 0.0;
      Tyy[t1][x][dimy-2] = -Tyy[t1][x][dimy-4];
      Tyy[t1][x][dimy-3] = 0.0;
      Txy[t1][x][dimy-2] = -Txy[t1][x][dimy-5];
      Txy[t1][x][dimy-3] = -Txy[t1][x][dimy-4];
    }
    #pragma omp for
    for(int y=0;y<dimy;y++){
      // boundary x=2
      Txx[t1][1][y] = -Txx[t1][3][y];
      Txx[t1][2][y] = 0.0;
      Tyy[t1][1][y] = -Tyy[t1][3][y];
      Tyy[t1][2][y] = 0.0;
      Txy[t1][0][y] = -Txy[t1][3][y];
      Txy[t1][1][y] = -Txy[t1][2][y];
      // boundary x=dimx+2
      Txx[t1][dimx-2][y] = -Txx[t1][dimx-4][y];
      Txx[t1][dimx-3][y] = 0.0;
      Tyy[t1][dimx-2][y] = -Tyy[t1][dimx-4][y];
      Tyy[t1][dimx-3][y] = 0.0;
      Txy[t1][dimx-2][y] = -Txy[t1][dimx-5][y];
      Txy[t1][dimx-3][y] = -Txy[t1][dimx-4][y];
    }

    // Compute velocities
    #pragma omp for
    for(int x=margin;x<dimx-margin;x++){
      for(int y=margin;y<dimy-margin;y++){
        U[t1][x][y]=(1.0F/24.0F)*beta*dt*(-27*Txx[t1][x][y] + Txx[t1][x - 1][y] + 27*Txx[t1][x + 1][y] - Txx[t1][x + 2][y] + 27*Txy[t1][x][y] + Txy[t1][x][y - 2] - 27*Txy[t1][x][y - 1] - Txy[t1][x][y + 1])/h + U[t][x][y];
        V[t1][x][y]=(1.0F/24.0F)*beta*dt*(27*Txy[t1][x][y] + Txy[t1][x - 2][y] - 27*Txy[t1][x - 1][y] - Txy[t1][x + 1][y] - 27*Tyy[t1][x][y] + Tyy[t1][x][y - 1] + 27*Tyy[t1][x][y + 1] - Tyy[t1][x][y + 2])/h + V[t][x][y];
      }
    }

    // update ghost cells for boundary conditions
    #pragma omp for
    for(int x=0;x<dimx;x++){
      // boundary y=2
      U[t1][x][1]=U[t1][x][3] - V[t1][x][1] - V[t1][x][2] + V[t1][x + 1][1] + V[t1][x + 1][2];
      V[t1][x][1]=(lambda*U[t1][x][2] - lambda*U[t1][x - 1][2] + lambda*V[t1][x][2] + 2*mu*V[t1][x][2])/(lambda + 2*mu);
      // boundary y=dimy+2
      U[t1][x][dimy - 2]=U[t1][x][dimy - 4] + V[t1][x][dimy - 4] + V[t1][x][dimy - 3] - V[t1][x + 1][dimy - 4] - V[t1][x + 1][dimy - 3];
      V[t1][x][dimy - 3]=(-lambda*U[t1][x][dimy - 3] + lambda*U[t1][x - 1][dimy - 3] + lambda*V[t1][x][dimy - 4] + 2*mu*V[t1][x][dimy - 4])/(lambda + 2*mu);
    }
    #pragma omp for
    for(int y=0;y<dimy;y++){
      // boundary x=2
      U[t1][1][y]=(lambda*(U[t1][2][y] + V[t1][2][y] - V[t1][2][y - 1]) + 2*mu*V[t1][2][y] - 2*mu*V[t1][2][y - 1])/lambda;
      V[t1][1][y]=-U[t1][1][y] + U[t1][1][y + 1] - U[t1][2][y] + U[t1][2][y + 1] + V[t1][3][y];
      // boundary x=dimy-3
      U[t1][dimx - 3][y]=(lambda*(U[t1][dimx - 4][y] - V[t1][dimx - 3][y] + V[t1][dimx - 3][y - 1]) - 2*mu*V[t1][dimx - 3][y] + 2*mu*V[t1][dimx - 3][y - 1])/lambda;
      V[t1][dimx - 2][y]=U[t1][dimx - 4][y] - U[t1][dimx - 4][y + 1] + U[t1][dimx - 3][y] - U[t1][dimx - 3][y + 1] + V[t1][dimx - 4][y];
    }
  } // end of time loop
  } // end of parallel section

  return 0;
}
