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
#include <cstdio>
#include <string>

int main(){

  const int _tp = 2;
  int t1 = 0;
  int t = 0;
  
  const float dt = 0.005;
const float mu = 0.5;
const float rho = 1.0;
const float dz = 0.02;
const float dx = 0.02;
const float tmax = 4.0;
const float dy = 0.02;
const float beta = 1.0;
const float lambda = 0.5;
const int dimx = 55;
const int ntsteps = 800;
const int dimy = 55;
const int dimz = 55;

  std::vector<float> _Txx_vec(2*dimx*dimy*dimz);
float (*Txx)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _Txx_vec.data();
std::vector<float> _Tyy_vec(2*dimx*dimy*dimz);
float (*Tyy)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _Tyy_vec.data();
std::vector<float> _Tzz_vec(2*dimx*dimy*dimz);
float (*Tzz)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _Tzz_vec.data();
std::vector<float> _Txy_vec(2*dimx*dimy*dimz);
float (*Txy)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _Txy_vec.data();
std::vector<float> _Tyz_vec(2*dimx*dimy*dimz);
float (*Tyz)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _Tyz_vec.data();
std::vector<float> _Txz_vec(2*dimx*dimy*dimz);
float (*Txz)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _Txz_vec.data();
std::vector<float> _U_vec(2*dimx*dimy*dimz);
float (*U)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _U_vec.data();
std::vector<float> _V_vec(2*dimx*dimy*dimz);
float (*V)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _V_vec.data();
std::vector<float> _W_vec(2*dimx*dimy*dimz);
float (*W)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) _W_vec.data();


#pragma omp parallel
  {
  #pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 2);
    		float z = dz*(k - 2);
    		Txx[0][i][j][k]=0;
    	}
    }
}
#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 2);
    		float z = dz*(k - 2);
    		Tyy[0][i][j][k]=0;
    	}
    }
}
#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 2);
    		float z = dz*(k - 2);
    		Tzz[0][i][j][k]=0;
    	}
    }
}
#pragma omp for
for(int i=2;i<dimx - 3;++i){
    for(int j=2;j<dimy - 3;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 1.5);
    		float y = dy*(j - 1.5);
    		float z = dz*(k - 2);
    		Txy[0][i][j][k]=0.0;
    	}
    }
}
#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 3;++j){
    	for(int k=2;k<dimz - 3;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 1.5);
    		float z = dz*(k - 1.5);
    		Tyz[0][i][j][k]=0.0;
    	}
    }
}
#pragma omp for
for(int i=2;i<dimx - 3;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 3;++k){
    		float x = dx*(i - 1.5);
    		float y = dy*(j - 2);
    		float z = dz*(k - 1.5);
    		Txz[0][i][j][k]=0.0;
    	}
    }
}
#pragma omp for
for(int i=2;i<dimx - 3;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 1.5);
    		float y = dy*(j - 2);
    		float z = dz*(k - 2);
    		U[0][i][j][k]=(sin(M_PI*y) - sin(M_PI*z))*cos(M_PI*x)*cos(0.0025*sqrt(2)*M_PI*sqrt(mu/rho));
    	}
    }
}
#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 3;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 1.5);
    		float z = dz*(k - 2);
    		V[0][i][j][k]=(-sin(M_PI*x) + sin(M_PI*z))*cos(M_PI*y)*cos(0.0025*sqrt(2)*M_PI*sqrt(mu/rho));
    	}
    }
}
#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 3;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 2);
    		float z = dz*(k - 1.5);
    		W[0][i][j][k]=(sin(M_PI*x) - sin(M_PI*y))*cos(M_PI*z)*cos(0.0025*sqrt(2)*M_PI*sqrt(mu/rho));
    	}
    }
}

    // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txx[0][2][y][z] = 0;
				Txx[0][1][y][z] = -Txx[0][3][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txx[0][dimx - 3][y][z] = 0;
				Txx[0][dimx - 2][y][z] = -Txx[0][dimx - 4][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Tyy[0][x][2][z] = 0;
				Tyy[0][x][1][z] = -Tyy[0][x][3][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Tyy[0][x][dimy - 3][z] = 0;
				Tyy[0][x][dimy - 2][z] = -Tyy[0][x][dimy - 4][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Tzz[0][x][y][2] = 0;
				Tzz[0][x][y][1] = -Tzz[0][x][y][3];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Tzz[0][x][y][dimz - 3] = 0;
				Tzz[0][x][y][dimz - 2] = -Tzz[0][x][y][dimz - 4];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txy[0][1][y][z] = -Txy[0][2][y][z];
				Txy[0][0][y][z] = -Txy[0][3][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txy[0][dimx - 3][y][z] = -Txy[0][dimx - 4][y][z];
				Txy[0][dimx - 2][y][z] = -Txy[0][dimx - 5][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Txy[0][x][1][z] = -Txy[0][x][2][z];
				Txy[0][x][0][z] = -Txy[0][x][3][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Txy[0][x][dimy - 3][z] = -Txy[0][x][dimy - 4][z];
				Txy[0][x][dimy - 2][z] = -Txy[0][x][dimy - 5][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Tyz[0][x][1][z] = -Tyz[0][x][2][z];
				Tyz[0][x][0][z] = -Tyz[0][x][3][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Tyz[0][x][dimy - 3][z] = -Tyz[0][x][dimy - 4][z];
				Tyz[0][x][dimy - 2][z] = -Tyz[0][x][dimy - 5][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Tyz[0][x][y][1] = -Tyz[0][x][y][2];
				Tyz[0][x][y][0] = -Tyz[0][x][y][3];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Tyz[0][x][y][dimz - 3] = -Tyz[0][x][y][dimz - 4];
				Tyz[0][x][y][dimz - 2] = -Tyz[0][x][y][dimz - 5];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txz[0][1][y][z] = -Txz[0][2][y][z];
				Txz[0][0][y][z] = -Txz[0][3][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txz[0][dimx - 3][y][z] = -Txz[0][dimx - 4][y][z];
				Txz[0][dimx - 2][y][z] = -Txz[0][dimx - 5][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Txz[0][x][y][1] = -Txz[0][x][y][2];
				Txz[0][x][y][0] = -Txz[0][x][y][3];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Txz[0][x][y][dimz - 3] = -Txz[0][x][y][dimz - 4];
				Txz[0][x][y][dimz - 2] = -Txz[0][x][y][dimz - 5];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        U[0][1][y][z]=(dx*dy*lambda*W[0][2][y][z] - dx*dy*lambda*W[0][2][y][z - 1] + dx*dz*lambda*V[0][2][y][z] - dx*dz*lambda*V[0][2][y - 1][z] + dy*dz*lambda*U[0][2][y][z] + 2*dy*dz*mu*U[0][2][y][z])/(dy*dz*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        U[0][dimx - 3][y][z]=(-dx*dy*lambda*W[0][dimx - 3][y][z] + dx*dy*lambda*W[0][dimx - 3][y][z - 1] - dx*dz*lambda*V[0][dimx - 3][y][z] + dx*dz*lambda*V[0][dimx - 3][y - 1][z] + dy*dz*lambda*U[0][dimx - 4][y][z] + 2*dy*dz*mu*U[0][dimx - 4][y][z])/(dy*dz*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        V[0][1][y][z]=(-dx*U[0][1][y][z] + dx*U[0][1][y + 1][z] + dx*U[0][2][y][z] - dx*U[0][2][y + 1][z] + dy*(2*V[0][2][y][z] - V[0][3][y][z]))/dy;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        V[0][dimx - 2][y][z]=(-dx*U[0][dimx - 4][y][z] + dx*U[0][dimx - 4][y + 1][z] + dx*U[0][dimx - 3][y][z] - dx*U[0][dimx - 3][y + 1][z] + dy*(-V[0][dimx - 4][y][z] + 2*V[0][dimx - 3][y][z]))/dy;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        W[0][1][y][z]=(-dx*U[0][1][y][z] + dx*U[0][1][y][z + 1] + dx*U[0][2][y][z] - dx*U[0][2][y][z + 1] + dz*(2*W[0][2][y][z] - W[0][3][y][z]))/dz;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        W[0][dimx - 2][y][z]=(-dx*U[0][dimx - 4][y][z] + dx*U[0][dimx - 4][y][z + 1] + dx*U[0][dimx - 3][y][z] - dx*U[0][dimx - 3][y][z + 1] + dz*(-W[0][dimx - 4][y][z] + 2*W[0][dimx - 3][y][z]))/dz;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        V[0][x][1][z]=(dx*dy*lambda*W[0][x][2][z] - dx*dy*lambda*W[0][x][2][z - 1] + dx*dz*lambda*V[0][x][2][z] + 2*dx*dz*mu*V[0][x][2][z] + dy*dz*lambda*U[0][x][2][z] - dy*dz*lambda*U[0][x - 1][2][z])/(dx*dz*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        V[0][x][dimy - 3][z]=(-dx*dy*lambda*W[0][x][dimy - 3][z] + dx*dy*lambda*W[0][x][dimy - 3][z - 1] + dx*dz*lambda*V[0][x][dimy - 4][z] + 2*dx*dz*mu*V[0][x][dimy - 4][z] - dy*dz*lambda*U[0][x][dimy - 3][z] + dy*dz*lambda*U[0][x - 1][dimy - 3][z])/(dx*dz*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        U[0][x][1][z]=(dx*(2*U[0][x][2][z] - U[0][x][3][z]) - dy*V[0][x][1][z] + dy*V[0][x][2][z] + dy*V[0][x + 1][1][z] - dy*V[0][x + 1][2][z])/dx;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        U[0][x][dimy - 2][z]=(dx*(-U[0][x][dimy - 4][z] + 2*U[0][x][dimy - 3][z]) - dy*V[0][x][dimy - 4][z] + dy*V[0][x][dimy - 3][z] + dy*V[0][x + 1][dimy - 4][z] - dy*V[0][x + 1][dimy - 3][z])/dx;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        W[0][x][1][z]=(-dy*V[0][x][1][z] + dy*V[0][x][1][z + 1] + dy*V[0][x][2][z] - dy*V[0][x][2][z + 1] + dz*(2*W[0][x][2][z] - W[0][x][3][z]))/dz;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        W[0][x][dimy - 2][z]=(-dy*V[0][x][dimy - 4][z] + dy*V[0][x][dimy - 4][z + 1] + dy*V[0][x][dimy - 3][z] - dy*V[0][x][dimy - 3][z + 1] + dz*(-W[0][x][dimy - 4][z] + 2*W[0][x][dimy - 3][z]))/dz;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        W[0][x][y][1]=(dx*dy*lambda*W[0][x][y][2] + 2*dx*dy*mu*W[0][x][y][2] + dx*dz*lambda*V[0][x][y][2] - dx*dz*lambda*V[0][x][y - 1][2] + dy*dz*lambda*U[0][x][y][2] - dy*dz*lambda*U[0][x - 1][y][2])/(dx*dy*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        W[0][x][y][dimz - 3]=(dx*dy*lambda*W[0][x][y][dimz - 4] + 2*dx*dy*mu*W[0][x][y][dimz - 4] - dx*dz*lambda*V[0][x][y][dimz - 3] + dx*dz*lambda*V[0][x][y - 1][dimz - 3] - dy*dz*lambda*U[0][x][y][dimz - 3] + dy*dz*lambda*U[0][x - 1][y][dimz - 3])/(dx*dy*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        U[0][x][y][1]=(dx*(2*U[0][x][y][2] - U[0][x][y][3]) - dz*W[0][x][y][1] + dz*W[0][x][y][2] + dz*W[0][x + 1][y][1] - dz*W[0][x + 1][y][2])/dx;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        U[0][x][y][dimz - 2]=(dx*(-U[0][x][y][dimz - 4] + 2*U[0][x][y][dimz - 3]) - dz*W[0][x][y][dimz - 4] + dz*W[0][x][y][dimz - 3] + dz*W[0][x + 1][y][dimz - 4] - dz*W[0][x + 1][y][dimz - 3])/dx;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        V[0][x][y][1]=(dy*(2*V[0][x][y][2] - V[0][x][y][3]) - dz*W[0][x][y][1] + dz*W[0][x][y][2] + dz*W[0][x][y + 1][1] - dz*W[0][x][y + 1][2])/dy;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        V[0][x][y][dimz - 2]=(dy*(-V[0][x][y][dimz - 4] + 2*V[0][x][y][dimz - 3]) - dz*W[0][x][y][dimz - 4] + dz*W[0][x][y][dimz - 3] + dz*W[0][x][y + 1][dimz - 4] - dz*W[0][x][y + 1][dimz - 3])/dy;
				
      }
  }


  for(int _ti=0;_ti<ntsteps;_ti++){
    
    // shared variables
    #pragma omp single
    {
      t = _ti % _tp; // array index of current time step
      t1 = (t+1) % _tp; // array index of the grid to be updated
    }

    #pragma omp for
for(int x=2;x<dimx - 2;++x){
    for(int y=2;y<dimy - 2;++y){
    	for(int z=2;z<dimz - 2;++z){
    		Txx[t1][x][y][z]=(1.0F/24.0F)*(dt*dx*dy*lambda*(27*W[t][x][y][z] + W[t][x][y][z - 2] - 27*W[t][x][y][z - 1] - W[t][x][y][z + 1]) + dt*dx*dz*lambda*(27*V[t][x][y][z] + V[t][x][y - 2][z] - 27*V[t][x][y - 1][z] - V[t][x][y + 1][z]) + dt*dy*dz*(27*lambda*U[t][x][y][z] + lambda*U[t][x - 2][y][z] - 27*lambda*U[t][x - 1][y][z] - lambda*U[t][x + 1][y][z] + 54*mu*U[t][x][y][z] + 2*mu*U[t][x - 2][y][z] - 54*mu*U[t][x - 1][y][z] - 2*mu*U[t][x + 1][y][z]) + 24*dx*dy*dz*Txx[t][x][y][z])/(dx*dy*dz);
			Tyy[t1][x][y][z]=(1.0F/24.0F)*(dt*dx*dy*lambda*(27*W[t][x][y][z] + W[t][x][y][z - 2] - 27*W[t][x][y][z - 1] - W[t][x][y][z + 1]) + dt*dx*dz*(27*lambda*V[t][x][y][z] + lambda*V[t][x][y - 2][z] - 27*lambda*V[t][x][y - 1][z] - lambda*V[t][x][y + 1][z] + 54*mu*V[t][x][y][z] + 2*mu*V[t][x][y - 2][z] - 54*mu*V[t][x][y - 1][z] - 2*mu*V[t][x][y + 1][z]) + dt*dy*dz*lambda*(27*U[t][x][y][z] + U[t][x - 2][y][z] - 27*U[t][x - 1][y][z] - U[t][x + 1][y][z]) + 24*dx*dy*dz*Tyy[t][x][y][z])/(dx*dy*dz);
			Tzz[t1][x][y][z]=(1.0F/24.0F)*(dt*dx*dy*(27*lambda*W[t][x][y][z] + lambda*W[t][x][y][z - 2] - 27*lambda*W[t][x][y][z - 1] - lambda*W[t][x][y][z + 1] + 54*mu*W[t][x][y][z] + 2*mu*W[t][x][y][z - 2] - 54*mu*W[t][x][y][z - 1] - 2*mu*W[t][x][y][z + 1]) + dt*dx*dz*lambda*(27*V[t][x][y][z] + V[t][x][y - 2][z] - 27*V[t][x][y - 1][z] - V[t][x][y + 1][z]) + dt*dy*dz*lambda*(27*U[t][x][y][z] + U[t][x - 2][y][z] - 27*U[t][x - 1][y][z] - U[t][x + 1][y][z]) + 24*dx*dy*dz*Tzz[t][x][y][z])/(dx*dy*dz);
			Txy[t1][x][y][z]=(1.0F/24.0F)*(dt*dx*mu*(-27*U[t][x][y][z] + U[t][x][y - 1][z] + 27*U[t][x][y + 1][z] - U[t][x][y + 2][z]) + dt*dy*mu*(-27*V[t][x][y][z] + V[t][x - 1][y][z] + 27*V[t][x + 1][y][z] - V[t][x + 2][y][z]) + 24*dx*dy*Txy[t][x][y][z])/(dx*dy);
			Tyz[t1][x][y][z]=(1.0F/24.0F)*(dt*dy*mu*(-27*V[t][x][y][z] + V[t][x][y][z - 1] + 27*V[t][x][y][z + 1] - V[t][x][y][z + 2]) + dt*dz*mu*(-27*W[t][x][y][z] + W[t][x][y - 1][z] + 27*W[t][x][y + 1][z] - W[t][x][y + 2][z]) + 24*dy*dz*Tyz[t][x][y][z])/(dy*dz);
			Txz[t1][x][y][z]=(1.0F/24.0F)*(dt*dx*mu*(-27*U[t][x][y][z] + U[t][x][y][z - 1] + 27*U[t][x][y][z + 1] - U[t][x][y][z + 2]) + dt*dz*mu*(-27*W[t][x][y][z] + W[t][x - 1][y][z] + 27*W[t][x + 1][y][z] - W[t][x + 2][y][z]) + 24*dx*dz*Txz[t][x][y][z])/(dx*dz);
			
    	}    	
    }
}

      // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txx[t1][2][y][z] = 0;
				Txx[t1][1][y][z] = -Txx[t1][3][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txx[t1][dimx - 3][y][z] = 0;
				Txx[t1][dimx - 2][y][z] = -Txx[t1][dimx - 4][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Tyy[t1][x][2][z] = 0;
				Tyy[t1][x][1][z] = -Tyy[t1][x][3][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Tyy[t1][x][dimy - 3][z] = 0;
				Tyy[t1][x][dimy - 2][z] = -Tyy[t1][x][dimy - 4][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Tzz[t1][x][y][2] = 0;
				Tzz[t1][x][y][1] = -Tzz[t1][x][y][3];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Tzz[t1][x][y][dimz - 3] = 0;
				Tzz[t1][x][y][dimz - 2] = -Tzz[t1][x][y][dimz - 4];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txy[t1][1][y][z] = -Txy[t1][2][y][z];
				Txy[t1][0][y][z] = -Txy[t1][3][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txy[t1][dimx - 3][y][z] = -Txy[t1][dimx - 4][y][z];
				Txy[t1][dimx - 2][y][z] = -Txy[t1][dimx - 5][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Txy[t1][x][1][z] = -Txy[t1][x][2][z];
				Txy[t1][x][0][z] = -Txy[t1][x][3][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Txy[t1][x][dimy - 3][z] = -Txy[t1][x][dimy - 4][z];
				Txy[t1][x][dimy - 2][z] = -Txy[t1][x][dimy - 5][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Tyz[t1][x][1][z] = -Tyz[t1][x][2][z];
				Tyz[t1][x][0][z] = -Tyz[t1][x][3][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        Tyz[t1][x][dimy - 3][z] = -Tyz[t1][x][dimy - 4][z];
				Tyz[t1][x][dimy - 2][z] = -Tyz[t1][x][dimy - 5][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Tyz[t1][x][y][1] = -Tyz[t1][x][y][2];
				Tyz[t1][x][y][0] = -Tyz[t1][x][y][3];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Tyz[t1][x][y][dimz - 3] = -Tyz[t1][x][y][dimz - 4];
				Tyz[t1][x][y][dimz - 2] = -Tyz[t1][x][y][dimz - 5];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txz[t1][1][y][z] = -Txz[t1][2][y][z];
				Txz[t1][0][y][z] = -Txz[t1][3][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
        Txz[t1][dimx - 3][y][z] = -Txz[t1][dimx - 4][y][z];
				Txz[t1][dimx - 2][y][z] = -Txz[t1][dimx - 5][y][z];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int z=0;z<dimz;++z){
        // nothing
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Txz[t1][x][y][1] = -Txz[t1][x][y][2];
				Txz[t1][x][y][0] = -Txz[t1][x][y][3];
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=0;x<dimx;++x){
      for(int y=0;y<dimy;++y){
        Txz[t1][x][y][dimz - 3] = -Txz[t1][x][y][dimz - 4];
				Txz[t1][x][y][dimz - 2] = -Txz[t1][x][y][dimz - 5];
				
      }
  }

    #pragma omp for
for(int x=2;x<dimx - 2;++x){
    for(int y=2;y<dimy - 2;++y){
    	for(int z=2;z<dimz - 2;++z){
    		U[t1][x][y][z]=(1.0F/24.0F)*(beta*dt*dx*dy*(27*Txz[t1][x][y][z] + Txz[t1][x][y][z - 2] - 27*Txz[t1][x][y][z - 1] - Txz[t1][x][y][z + 1]) + beta*dt*dx*dz*(27*Txy[t1][x][y][z] + Txy[t1][x][y - 2][z] - 27*Txy[t1][x][y - 1][z] - Txy[t1][x][y + 1][z]) + beta*dt*dy*dz*(-27*Txx[t1][x][y][z] + Txx[t1][x - 1][y][z] + 27*Txx[t1][x + 1][y][z] - Txx[t1][x + 2][y][z]) + 24*dx*dy*dz*U[t][x][y][z])/(dx*dy*dz);
			V[t1][x][y][z]=(1.0F/24.0F)*(beta*dt*dx*dy*(27*Tyz[t1][x][y][z] + Tyz[t1][x][y][z - 2] - 27*Tyz[t1][x][y][z - 1] - Tyz[t1][x][y][z + 1]) + beta*dt*dx*dz*(-27*Tyy[t1][x][y][z] + Tyy[t1][x][y - 1][z] + 27*Tyy[t1][x][y + 1][z] - Tyy[t1][x][y + 2][z]) + beta*dt*dy*dz*(27*Txy[t1][x][y][z] + Txy[t1][x - 2][y][z] - 27*Txy[t1][x - 1][y][z] - Txy[t1][x + 1][y][z]) + 24*dx*dy*dz*V[t][x][y][z])/(dx*dy*dz);
			W[t1][x][y][z]=(1.0F/24.0F)*(beta*dt*dx*dy*(-27*Tzz[t1][x][y][z] + Tzz[t1][x][y][z - 1] + 27*Tzz[t1][x][y][z + 1] - Tzz[t1][x][y][z + 2]) + beta*dt*dx*dz*(27*Tyz[t1][x][y][z] + Tyz[t1][x][y - 2][z] - 27*Tyz[t1][x][y - 1][z] - Tyz[t1][x][y + 1][z]) + beta*dt*dy*dz*(27*Txz[t1][x][y][z] + Txz[t1][x - 2][y][z] - 27*Txz[t1][x - 1][y][z] - Txz[t1][x + 1][y][z]) + 24*dx*dy*dz*W[t][x][y][z])/(dx*dy*dz);
			
    	}    	
    }
}

      // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        U[t1][1][y][z]=(dx*dy*lambda*W[t1][2][y][z] - dx*dy*lambda*W[t1][2][y][z - 1] + dx*dz*lambda*V[t1][2][y][z] - dx*dz*lambda*V[t1][2][y - 1][z] + dy*dz*lambda*U[t1][2][y][z] + 2*dy*dz*mu*U[t1][2][y][z])/(dy*dz*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        U[t1][dimx - 3][y][z]=(-dx*dy*lambda*W[t1][dimx - 3][y][z] + dx*dy*lambda*W[t1][dimx - 3][y][z - 1] - dx*dz*lambda*V[t1][dimx - 3][y][z] + dx*dz*lambda*V[t1][dimx - 3][y - 1][z] + dy*dz*lambda*U[t1][dimx - 4][y][z] + 2*dy*dz*mu*U[t1][dimx - 4][y][z])/(dy*dz*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        V[t1][1][y][z]=(-dx*U[t1][1][y][z] + dx*U[t1][1][y + 1][z] + dx*U[t1][2][y][z] - dx*U[t1][2][y + 1][z] + dy*(2*V[t1][2][y][z] - V[t1][3][y][z]))/dy;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        V[t1][dimx - 2][y][z]=(-dx*U[t1][dimx - 4][y][z] + dx*U[t1][dimx - 4][y + 1][z] + dx*U[t1][dimx - 3][y][z] - dx*U[t1][dimx - 3][y + 1][z] + dy*(-V[t1][dimx - 4][y][z] + 2*V[t1][dimx - 3][y][z]))/dy;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        W[t1][1][y][z]=(-dx*U[t1][1][y][z] + dx*U[t1][1][y][z + 1] + dx*U[t1][2][y][z] - dx*U[t1][2][y][z + 1] + dz*(2*W[t1][2][y][z] - W[t1][3][y][z]))/dz;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int y=1;y<dimy - 1;++y){
      for(int z=1;z<dimz - 1;++z){
        W[t1][dimx - 2][y][z]=(-dx*U[t1][dimx - 4][y][z] + dx*U[t1][dimx - 4][y][z + 1] + dx*U[t1][dimx - 3][y][z] - dx*U[t1][dimx - 3][y][z + 1] + dz*(-W[t1][dimx - 4][y][z] + 2*W[t1][dimx - 3][y][z]))/dz;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        V[t1][x][1][z]=(dx*dy*lambda*W[t1][x][2][z] - dx*dy*lambda*W[t1][x][2][z - 1] + dx*dz*lambda*V[t1][x][2][z] + 2*dx*dz*mu*V[t1][x][2][z] + dy*dz*lambda*U[t1][x][2][z] - dy*dz*lambda*U[t1][x - 1][2][z])/(dx*dz*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        V[t1][x][dimy - 3][z]=(-dx*dy*lambda*W[t1][x][dimy - 3][z] + dx*dy*lambda*W[t1][x][dimy - 3][z - 1] + dx*dz*lambda*V[t1][x][dimy - 4][z] + 2*dx*dz*mu*V[t1][x][dimy - 4][z] - dy*dz*lambda*U[t1][x][dimy - 3][z] + dy*dz*lambda*U[t1][x - 1][dimy - 3][z])/(dx*dz*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        U[t1][x][1][z]=(dx*(2*U[t1][x][2][z] - U[t1][x][3][z]) - dy*V[t1][x][1][z] + dy*V[t1][x][2][z] + dy*V[t1][x + 1][1][z] - dy*V[t1][x + 1][2][z])/dx;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        U[t1][x][dimy - 2][z]=(dx*(-U[t1][x][dimy - 4][z] + 2*U[t1][x][dimy - 3][z]) - dy*V[t1][x][dimy - 4][z] + dy*V[t1][x][dimy - 3][z] + dy*V[t1][x + 1][dimy - 4][z] - dy*V[t1][x + 1][dimy - 3][z])/dx;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        W[t1][x][1][z]=(-dy*V[t1][x][1][z] + dy*V[t1][x][1][z + 1] + dy*V[t1][x][2][z] - dy*V[t1][x][2][z + 1] + dz*(2*W[t1][x][2][z] - W[t1][x][3][z]))/dz;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int z=1;z<dimz - 1;++z){
        W[t1][x][dimy - 2][z]=(-dy*V[t1][x][dimy - 4][z] + dy*V[t1][x][dimy - 4][z + 1] + dy*V[t1][x][dimy - 3][z] - dy*V[t1][x][dimy - 3][z + 1] + dz*(-W[t1][x][dimy - 4][z] + 2*W[t1][x][dimy - 3][z]))/dz;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        W[t1][x][y][1]=(dx*dy*lambda*W[t1][x][y][2] + 2*dx*dy*mu*W[t1][x][y][2] + dx*dz*lambda*V[t1][x][y][2] - dx*dz*lambda*V[t1][x][y - 1][2] + dy*dz*lambda*U[t1][x][y][2] - dy*dz*lambda*U[t1][x - 1][y][2])/(dx*dy*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        W[t1][x][y][dimz - 3]=(dx*dy*lambda*W[t1][x][y][dimz - 4] + 2*dx*dy*mu*W[t1][x][y][dimz - 4] - dx*dz*lambda*V[t1][x][y][dimz - 3] + dx*dz*lambda*V[t1][x][y - 1][dimz - 3] - dy*dz*lambda*U[t1][x][y][dimz - 3] + dy*dz*lambda*U[t1][x - 1][y][dimz - 3])/(dx*dy*(lambda + 2*mu));
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        U[t1][x][y][1]=(dx*(2*U[t1][x][y][2] - U[t1][x][y][3]) - dz*W[t1][x][y][1] + dz*W[t1][x][y][2] + dz*W[t1][x + 1][y][1] - dz*W[t1][x + 1][y][2])/dx;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        U[t1][x][y][dimz - 2]=(dx*(-U[t1][x][y][dimz - 4] + 2*U[t1][x][y][dimz - 3]) - dz*W[t1][x][y][dimz - 4] + dz*W[t1][x][y][dimz - 3] + dz*W[t1][x + 1][y][dimz - 4] - dz*W[t1][x + 1][y][dimz - 3])/dx;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        V[t1][x][y][1]=(dy*(2*V[t1][x][y][2] - V[t1][x][y][3]) - dz*W[t1][x][y][1] + dz*W[t1][x][y][2] + dz*W[t1][x][y + 1][1] - dz*W[t1][x][y + 1][2])/dy;
				
      }
  }
  // update ghost cells for boundary conditions
  #pragma omp for
  for(int x=1;x<dimx - 1;++x){
      for(int y=1;y<dimy - 1;++y){
        V[t1][x][y][dimz - 2]=(dy*(-V[t1][x][y][dimz - 4] + 2*V[t1][x][y][dimz - 3]) - dz*W[t1][x][y][dimz - 4] + dz*W[t1][x][y][dimz - 3] + dz*W[t1][x][y + 1][dimz - 4] - dz*W[t1][x][y + 1][dimz - 3])/dy;
				
      }
  }


    

  } // end of time loop
  } // end of parallel section

  float Txx_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 2);
    		float z = dz*(k - 2);
    		Txx_l2+=pow(sqrt(2)*sqrt(mu*rho)*(sin(M_PI*y) - sin(M_PI*z))*sin(M_PI*x)*sin(4.0*sqrt(2)*M_PI*sqrt(mu/rho)) + Txx[0][i][j][k], 2.0);
    	}
    }
}
printf("Txx_l2 = %.10f\n", Txx_l2);
		float Tyy_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 2);
    		float z = dz*(k - 2);
    		Tyy_l2+=pow(sqrt(2)*sqrt(mu*rho)*(-sin(M_PI*x) + sin(M_PI*z))*sin(M_PI*y)*sin(4.0*sqrt(2)*M_PI*sqrt(mu/rho)) + Tyy[0][i][j][k], 2.0);
    	}
    }
}
printf("Tyy_l2 = %.10f\n", Tyy_l2);
		float Tzz_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 2);
    		float z = dz*(k - 2);
    		Tzz_l2+=pow(sqrt(2)*sqrt(mu*rho)*(sin(M_PI*x) - sin(M_PI*y))*sin(M_PI*z)*sin(4.0*sqrt(2)*M_PI*sqrt(mu/rho)) + Tzz[0][i][j][k], 2.0);
    	}
    }
}
printf("Tzz_l2 = %.10f\n", Tzz_l2);
		float Txy_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 3;++i){
    for(int j=2;j<dimy - 3;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 1.5);
    		float y = dy*(j - 1.5);
    		float z = dz*(k - 2);
    		Txy_l2+=pow(Txy[0][i][j][k], 2.0);
    	}
    }
}
printf("Txy_l2 = %.10f\n", Txy_l2);
		float Tyz_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 3;++j){
    	for(int k=2;k<dimz - 3;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 1.5);
    		float z = dz*(k - 1.5);
    		Tyz_l2+=pow(Tyz[0][i][j][k], 2.0);
    	}
    }
}
printf("Tyz_l2 = %.10f\n", Tyz_l2);
		float Txz_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 3;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 3;++k){
    		float x = dx*(i - 1.5);
    		float y = dy*(j - 2);
    		float z = dz*(k - 1.5);
    		Txz_l2+=pow(Txz[0][i][j][k], 2.0);
    	}
    }
}
printf("Txz_l2 = %.10f\n", Txz_l2);
		float U_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 3;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 1.5);
    		float y = dy*(j - 2);
    		float z = dz*(k - 2);
    		U_l2+=pow(-(sin(M_PI*y) - sin(M_PI*z))*cos(M_PI*x)*cos(4.0025*sqrt(2)*M_PI*sqrt(mu/rho)) + U[0][i][j][k], 2.0);
    	}
    }
}
printf("U_l2 = %.10f\n", U_l2);
		float V_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 3;++j){
    	for(int k=2;k<dimz - 2;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 1.5);
    		float z = dz*(k - 2);
    		V_l2+=pow(-(-sin(M_PI*x) + sin(M_PI*z))*cos(M_PI*y)*cos(4.0025*sqrt(2)*M_PI*sqrt(mu/rho)) + V[0][i][j][k], 2.0);
    	}
    }
}
printf("V_l2 = %.10f\n", V_l2);
		float W_l2 = 0.0;
		#pragma omp for
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
    	for(int k=2;k<dimz - 3;++k){
    		float x = dx*(i - 2);
    		float y = dy*(j - 2);
    		float z = dz*(k - 1.5);
    		W_l2+=pow(-(sin(M_PI*x) - sin(M_PI*y))*cos(M_PI*z)*cos(4.0025*sqrt(2)*M_PI*sqrt(mu/rho)) + W[0][i][j][k], 2.0);
    	}
    }
}
printf("W_l2 = %.10f\n", W_l2);
		

  return 0;
}
