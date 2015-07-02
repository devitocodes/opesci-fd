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


#include "opesciIO.h"
#include "opesciHandy.h"

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
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	float x = h*(i - 2);
    	float y = h*(j - 2);
    	Txx[0][i][j]=0
    }
}
#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	float x = h*(i - 2);
    	float y = h*(j - 2);
    	Tyy[0][i][j]=0
    }
}
#pragma omp for
for(int i=2;i<dimx - 3;++i){
    for(int j=2;j<dimy - 3;++j){
    	float x = h*(i - 1.5);
    	float y = h*(j - 1.5);
    	Txy[0][i][j]=0.0
    }
}
#pragma omp for
for(int i=2;i<dimx - 3;++i){
    for(int j=2;j<dimy - 2;++j){
    	float x = h*(i - 1.5);
    	float y = h*(j - 2);
    	U[0][i][j]=1.4142135623731*M_PI*Vs*sin(M_PI*y)*cos(M_PI*x)
    }
}
#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 3;++j){
    	float x = h*(i - 2);
    	float y = h*(j - 1.5);
    	V[0][i][j]=-1.4142135623731*M_PI*Vs*sin(M_PI*x)*cos(M_PI*y)
    }
}


  } // end of parallel section

  return 0;
}
