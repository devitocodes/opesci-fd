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

#ifdef _MSC_VER
#define M_PI 3.14159265358979323846
#endif
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

const int dim1 = 105;
const int dim2 = 105;
const int dim3 = 105;
const float dx1 = 0.01;
const float dx2 = 0.01;
const float dx3 = 0.01;
int t0 = 0;
int t1 = 0;
const int tp = 2;
const float dt = 0.002;
const int margin = 2;
const int ntsteps = 500;
const int vec_size = 2*dim1*dim2*dim3;
const float mu = 0.25;
const float beta = 1.0;
const float rho = 1.0;
const float lambda = 0.5;

float *_Txx_vec;
#ifdef _MSC_VER
_Txx_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_Txx_vec), 4096, 2315250*sizeof(float));
#endif
float (*Txx)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _Txx_vec;
float *_Tyy_vec;
#ifdef _MSC_VER
_Tyy_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_Tyy_vec), 4096, 2315250*sizeof(float));
#endif
float (*Tyy)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _Tyy_vec;
float *_Tzz_vec;
#ifdef _MSC_VER
_Tzz_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_Tzz_vec), 4096, 2315250*sizeof(float));
#endif
float (*Tzz)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _Tzz_vec;
float *_Txy_vec;
#ifdef _MSC_VER
_Txy_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_Txy_vec), 4096, 2315250*sizeof(float));
#endif
float (*Txy)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _Txy_vec;
float *_Tyz_vec;
#ifdef _MSC_VER
_Tyz_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_Tyz_vec), 4096, 2315250*sizeof(float));
#endif
float (*Tyz)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _Tyz_vec;
float *_Txz_vec;
#ifdef _MSC_VER
_Txz_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_Txz_vec), 4096, 2315250*sizeof(float));
#endif
float (*Txz)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _Txz_vec;
float *_U_vec;
#ifdef _MSC_VER
_U_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_U_vec), 4096, 2315250*sizeof(float));
#endif
float (*U)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _U_vec;
float *_V_vec;
#ifdef _MSC_VER
_V_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_V_vec), 4096, 2315250*sizeof(float));
#endif
float (*V)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _V_vec;
float *_W_vec;
#ifdef _MSC_VER
_W_vec = (float*) _aligned_malloc(2315250*sizeof(float), 4096);
#else
posix_memalign((void **)(&_W_vec), 4096, 2315250*sizeof(float));
#endif
float (*W)[dim1][dim2][dim3]= (float (*)[dim1][dim2][dim3]) _W_vec;


#pragma omp parallel
{
#pragma omp for
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Txx[0][_x][_y][_z]=0;
}
}
}
#pragma omp for
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Tyy[0][_x][_y][_z]=0;
}
}
}
#pragma omp for
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Tzz[0][_x][_y][_z]=0;
}
}
}
#pragma omp for
for(int _x=2;_x<dim1 - 3;++_x){
float x= dx1*(_x - 1.5);
for(int _y=2;_y<dim2 - 3;++_y){
float y= dx2*(_y - 1.5);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Txy[0][_x][_y][_z]=0.0;
}
}
}
#pragma omp for
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 3;++_y){
float y= dx2*(_y - 1.5);
for(int _z=2;_z<dim3 - 3;++_z){
float z= dx3*(_z - 1.5);
Tyz[0][_x][_y][_z]=0.0;
}
}
}
#pragma omp for
for(int _x=2;_x<dim1 - 3;++_x){
float x= dx1*(_x - 1.5);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 3;++_z){
float z= dx3*(_z - 1.5);
Txz[0][_x][_y][_z]=0.0;
}
}
}
#pragma omp for
for(int _x=2;_x<dim1 - 3;++_x){
float x= dx1*(_x - 1.5);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
U[0][_x][_y][_z]=(sin(M_PI*y) - sin(M_PI*z))*cos(M_PI*x)*cos(0.001*sqrt(2)*M_PI*sqrt(beta*mu));
}
}
}
#pragma omp for
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 3;++_y){
float y= dx2*(_y - 1.5);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
V[0][_x][_y][_z]=(-sin(M_PI*x) + sin(M_PI*z))*cos(M_PI*y)*cos(0.001*sqrt(2)*M_PI*sqrt(beta*mu));
}
}
}
#pragma omp for
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 3;++_z){
float z= dx3*(_z - 1.5);
W[0][_x][_y][_z]=(sin(M_PI*x) - sin(M_PI*y))*cos(M_PI*z)*cos(0.001*sqrt(2)*M_PI*sqrt(beta*mu));
}
}
}

#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txx[0][2][y][z] = 0;
Txx[0][1][y][z] = -Txx[0][3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txx[0][dim1 - 3][y][z] = 0;
Txx[0][dim1 - 2][y][z] = -Txx[0][dim1 - 4][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyy[0][x][2][z] = 0;
Tyy[0][x][1][z] = -Tyy[0][x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyy[0][x][dim2 - 3][z] = 0;
Tyy[0][x][dim2 - 2][z] = -Tyy[0][x][dim2 - 4][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tzz[0][x][y][2] = 0;
Tzz[0][x][y][1] = -Tzz[0][x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tzz[0][x][y][dim3 - 3] = 0;
Tzz[0][x][y][dim3 - 2] = -Tzz[0][x][y][dim3 - 4];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[0][1][y][z] = -Txy[0][2][y][z];
Txy[0][0][y][z] = -Txy[0][3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[0][dim1 - 3][y][z] = -Txy[0][dim1 - 4][y][z];
Txy[0][dim1 - 2][y][z] = -Txy[0][dim1 - 5][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[0][x][1][z] = -Txy[0][x][2][z];
Txy[0][x][0][z] = -Txy[0][x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[0][x][dim2 - 3][z] = -Txy[0][x][dim2 - 4][z];
Txy[0][x][dim2 - 2][z] = -Txy[0][x][dim2 - 5][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyz[0][x][1][z] = -Tyz[0][x][2][z];
Tyz[0][x][0][z] = -Tyz[0][x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyz[0][x][dim2 - 3][z] = -Tyz[0][x][dim2 - 4][z];
Tyz[0][x][dim2 - 2][z] = -Tyz[0][x][dim2 - 5][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tyz[0][x][y][1] = -Tyz[0][x][y][2];
Tyz[0][x][y][0] = -Tyz[0][x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tyz[0][x][y][dim3 - 3] = -Tyz[0][x][y][dim3 - 4];
Tyz[0][x][y][dim3 - 2] = -Tyz[0][x][y][dim3 - 5];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txz[0][1][y][z] = -Txz[0][2][y][z];
Txz[0][0][y][z] = -Txz[0][3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txz[0][dim1 - 3][y][z] = -Txz[0][dim1 - 4][y][z];
Txz[0][dim1 - 2][y][z] = -Txz[0][dim1 - 5][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Txz[0][x][y][1] = -Txz[0][x][y][2];
Txz[0][x][y][0] = -Txz[0][x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Txz[0][x][y][dim3 - 3] = -Txz[0][x][y][dim3 - 4];
Txz[0][x][y][dim3 - 2] = -Txz[0][x][y][dim3 - 5];
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[0][1][y][z] = (lambda*dx1*dx2*W[0][2][y][z] - lambda*dx1*dx2*W[0][2][y][z - 1] + lambda*dx1*dx3*V[0][2][y][z] - lambda*dx1*dx3*V[0][2][y - 1][z] + lambda*dx2*dx3*U[0][2][y][z] + 2*mu*dx2*dx3*U[0][2][y][z])/(dx2*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[0][dim1 - 3][y][z] = (-lambda*dx1*dx2*W[0][dim1 - 3][y][z] + lambda*dx1*dx2*W[0][dim1 - 3][y][z - 1] - lambda*dx1*dx3*V[0][dim1 - 3][y][z] + lambda*dx1*dx3*V[0][dim1 - 3][y - 1][z] + lambda*dx2*dx3*U[0][dim1 - 4][y][z] + 2*mu*dx2*dx3*U[0][dim1 - 4][y][z])/(dx2*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[0][1][y][z] = (-dx1*U[0][1][y][z] + dx1*U[0][1][y + 1][z] + dx1*U[0][2][y][z] - dx1*U[0][2][y + 1][z] + dx2*(2*V[0][2][y][z] - V[0][3][y][z]))/dx2;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[0][dim1 - 2][y][z] = (-dx1*U[0][dim1 - 4][y][z] + dx1*U[0][dim1 - 4][y + 1][z] + dx1*U[0][dim1 - 3][y][z] - dx1*U[0][dim1 - 3][y + 1][z] + dx2*(-V[0][dim1 - 4][y][z] + 2*V[0][dim1 - 3][y][z]))/dx2;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[0][1][y][z] = (-dx1*U[0][1][y][z] + dx1*U[0][1][y][z + 1] + dx1*U[0][2][y][z] - dx1*U[0][2][y][z + 1] + dx3*(2*W[0][2][y][z] - W[0][3][y][z]))/dx3;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[0][dim1 - 2][y][z] = (-dx1*U[0][dim1 - 4][y][z] + dx1*U[0][dim1 - 4][y][z + 1] + dx1*U[0][dim1 - 3][y][z] - dx1*U[0][dim1 - 3][y][z + 1] + dx3*(-W[0][dim1 - 4][y][z] + 2*W[0][dim1 - 3][y][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[0][x][1][z] = (lambda*dx1*dx2*W[0][x][2][z] - lambda*dx1*dx2*W[0][x][2][z - 1] + lambda*dx1*dx3*V[0][x][2][z] + lambda*dx2*dx3*U[0][x][2][z] - lambda*dx2*dx3*U[0][x - 1][2][z] + 2*mu*dx1*dx3*V[0][x][2][z])/(dx1*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[0][x][dim2 - 3][z] = (-lambda*dx1*dx2*W[0][x][dim2 - 3][z] + lambda*dx1*dx2*W[0][x][dim2 - 3][z - 1] + lambda*dx1*dx3*V[0][x][dim2 - 4][z] - lambda*dx2*dx3*U[0][x][dim2 - 3][z] + lambda*dx2*dx3*U[0][x - 1][dim2 - 3][z] + 2*mu*dx1*dx3*V[0][x][dim2 - 4][z])/(dx1*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[0][x][1][z] = (dx1*(2*U[0][x][2][z] - U[0][x][3][z]) - dx2*V[0][x][1][z] + dx2*V[0][x][2][z] + dx2*V[0][x + 1][1][z] - dx2*V[0][x + 1][2][z])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[0][x][dim2 - 2][z] = (dx1*(-U[0][x][dim2 - 4][z] + 2*U[0][x][dim2 - 3][z]) - dx2*V[0][x][dim2 - 4][z] + dx2*V[0][x][dim2 - 3][z] + dx2*V[0][x + 1][dim2 - 4][z] - dx2*V[0][x + 1][dim2 - 3][z])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[0][x][1][z] = (-dx2*V[0][x][1][z] + dx2*V[0][x][1][z + 1] + dx2*V[0][x][2][z] - dx2*V[0][x][2][z + 1] + dx3*(2*W[0][x][2][z] - W[0][x][3][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[0][x][dim2 - 2][z] = (-dx2*V[0][x][dim2 - 4][z] + dx2*V[0][x][dim2 - 4][z + 1] + dx2*V[0][x][dim2 - 3][z] - dx2*V[0][x][dim2 - 3][z + 1] + dx3*(-W[0][x][dim2 - 4][z] + 2*W[0][x][dim2 - 3][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
W[0][x][y][1] = (lambda*dx1*dx2*W[0][x][y][2] + lambda*dx1*dx3*V[0][x][y][2] - lambda*dx1*dx3*V[0][x][y - 1][2] + lambda*dx2*dx3*U[0][x][y][2] - lambda*dx2*dx3*U[0][x - 1][y][2] + 2*mu*dx1*dx2*W[0][x][y][2])/(dx1*dx2*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
W[0][x][y][dim3 - 3] = (lambda*dx1*dx2*W[0][x][y][dim3 - 4] - lambda*dx1*dx3*V[0][x][y][dim3 - 3] + lambda*dx1*dx3*V[0][x][y - 1][dim3 - 3] - lambda*dx2*dx3*U[0][x][y][dim3 - 3] + lambda*dx2*dx3*U[0][x - 1][y][dim3 - 3] + 2*mu*dx1*dx2*W[0][x][y][dim3 - 4])/(dx1*dx2*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
U[0][x][y][1] = (dx1*(2*U[0][x][y][2] - U[0][x][y][3]) - dx3*W[0][x][y][1] + dx3*W[0][x][y][2] + dx3*W[0][x + 1][y][1] - dx3*W[0][x + 1][y][2])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
U[0][x][y][dim3 - 2] = (dx1*(-U[0][x][y][dim3 - 4] + 2*U[0][x][y][dim3 - 3]) - dx3*W[0][x][y][dim3 - 4] + dx3*W[0][x][y][dim3 - 3] + dx3*W[0][x + 1][y][dim3 - 4] - dx3*W[0][x + 1][y][dim3 - 3])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
V[0][x][y][1] = (dx2*(2*V[0][x][y][2] - V[0][x][y][3]) - dx3*W[0][x][y][1] + dx3*W[0][x][y][2] + dx3*W[0][x][y + 1][1] - dx3*W[0][x][y + 1][2])/dx2;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
V[0][x][y][dim3 - 2] = (dx2*(-V[0][x][y][dim3 - 4] + 2*V[0][x][y][dim3 - 3]) - dx3*W[0][x][y][dim3 - 4] + dx3*W[0][x][y][dim3 - 3] + dx3*W[0][x][y + 1][dim3 - 4] - dx3*W[0][x][y + 1][dim3 - 3])/dx2;
}
}


for(int _ti=0;_ti<ntsteps;_ti++){

    #pragma omp single
    {
    	t0 = (_ti)%(tp);
t1 = (t0 + 1)%(tp);

    }

#pragma omp for
for(int x=2;x<dim1 - 2;++x){
for(int y=2;y<dim2 - 2;++y){
#pragma ivdep
for(int z=2;z<dim3 - 2;++z){
Txx[t1][x][y][z]=(1.0F/24.0F)*(lambda*dt*dx1*dx2*(27*W[t0][x][y][z] + W[t0][x][y][z - 2] - 27*W[t0][x][y][z - 1] - W[t0][x][y][z + 1]) + lambda*dt*dx1*dx3*(27*V[t0][x][y][z] + V[t0][x][y - 2][z] - 27*V[t0][x][y - 1][z] - V[t0][x][y + 1][z]) + dt*dx2*dx3*(27*lambda*U[t0][x][y][z] + lambda*U[t0][x - 2][y][z] - 27*lambda*U[t0][x - 1][y][z] - lambda*U[t0][x + 1][y][z] + 54*mu*U[t0][x][y][z] + 2*mu*U[t0][x - 2][y][z] - 54*mu*U[t0][x - 1][y][z] - 2*mu*U[t0][x + 1][y][z]) + 24*dx1*dx2*dx3*Txx[t0][x][y][z])/(dx1*dx2*dx3);
Tyy[t1][x][y][z]=(1.0F/24.0F)*(lambda*dt*dx1*dx2*(27*W[t0][x][y][z] + W[t0][x][y][z - 2] - 27*W[t0][x][y][z - 1] - W[t0][x][y][z + 1]) + lambda*dt*dx2*dx3*(27*U[t0][x][y][z] + U[t0][x - 2][y][z] - 27*U[t0][x - 1][y][z] - U[t0][x + 1][y][z]) + dt*dx1*dx3*(27*lambda*V[t0][x][y][z] + lambda*V[t0][x][y - 2][z] - 27*lambda*V[t0][x][y - 1][z] - lambda*V[t0][x][y + 1][z] + 54*mu*V[t0][x][y][z] + 2*mu*V[t0][x][y - 2][z] - 54*mu*V[t0][x][y - 1][z] - 2*mu*V[t0][x][y + 1][z]) + 24*dx1*dx2*dx3*Tyy[t0][x][y][z])/(dx1*dx2*dx3);
Tzz[t1][x][y][z]=(1.0F/24.0F)*(lambda*dt*dx1*dx3*(27*V[t0][x][y][z] + V[t0][x][y - 2][z] - 27*V[t0][x][y - 1][z] - V[t0][x][y + 1][z]) + lambda*dt*dx2*dx3*(27*U[t0][x][y][z] + U[t0][x - 2][y][z] - 27*U[t0][x - 1][y][z] - U[t0][x + 1][y][z]) + dt*dx1*dx2*(27*lambda*W[t0][x][y][z] + lambda*W[t0][x][y][z - 2] - 27*lambda*W[t0][x][y][z - 1] - lambda*W[t0][x][y][z + 1] + 54*mu*W[t0][x][y][z] + 2*mu*W[t0][x][y][z - 2] - 54*mu*W[t0][x][y][z - 1] - 2*mu*W[t0][x][y][z + 1]) + 24*dx1*dx2*dx3*Tzz[t0][x][y][z])/(dx1*dx2*dx3);
Txy[t1][x][y][z]=(1.0F/24.0F)*(mu*dt*dx1*(-27*U[t0][x][y][z] + U[t0][x][y - 1][z] + 27*U[t0][x][y + 1][z] - U[t0][x][y + 2][z]) + mu*dt*dx2*(-27*V[t0][x][y][z] + V[t0][x - 1][y][z] + 27*V[t0][x + 1][y][z] - V[t0][x + 2][y][z]) + 24*dx1*dx2*Txy[t0][x][y][z])/(dx1*dx2);
Tyz[t1][x][y][z]=(1.0F/24.0F)*(mu*dt*dx2*(-27*V[t0][x][y][z] + V[t0][x][y][z - 1] + 27*V[t0][x][y][z + 1] - V[t0][x][y][z + 2]) + mu*dt*dx3*(-27*W[t0][x][y][z] + W[t0][x][y - 1][z] + 27*W[t0][x][y + 1][z] - W[t0][x][y + 2][z]) + 24*dx2*dx3*Tyz[t0][x][y][z])/(dx2*dx3);
Txz[t1][x][y][z]=(1.0F/24.0F)*(mu*dt*dx1*(-27*U[t0][x][y][z] + U[t0][x][y][z - 1] + 27*U[t0][x][y][z + 1] - U[t0][x][y][z + 2]) + mu*dt*dx3*(-27*W[t0][x][y][z] + W[t0][x - 1][y][z] + 27*W[t0][x + 1][y][z] - W[t0][x + 2][y][z]) + 24*dx1*dx3*Txz[t0][x][y][z])/(dx1*dx3);
}
}
}

#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txx[t1][2][y][z] = 0;
Txx[t1][1][y][z] = -Txx[t1][3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txx[t1][dim1 - 3][y][z] = 0;
Txx[t1][dim1 - 2][y][z] = -Txx[t1][dim1 - 4][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyy[t1][x][2][z] = 0;
Tyy[t1][x][1][z] = -Tyy[t1][x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyy[t1][x][dim2 - 3][z] = 0;
Tyy[t1][x][dim2 - 2][z] = -Tyy[t1][x][dim2 - 4][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tzz[t1][x][y][2] = 0;
Tzz[t1][x][y][1] = -Tzz[t1][x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tzz[t1][x][y][dim3 - 3] = 0;
Tzz[t1][x][y][dim3 - 2] = -Tzz[t1][x][y][dim3 - 4];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[t1][1][y][z] = -Txy[t1][2][y][z];
Txy[t1][0][y][z] = -Txy[t1][3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[t1][dim1 - 3][y][z] = -Txy[t1][dim1 - 4][y][z];
Txy[t1][dim1 - 2][y][z] = -Txy[t1][dim1 - 5][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[t1][x][1][z] = -Txy[t1][x][2][z];
Txy[t1][x][0][z] = -Txy[t1][x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[t1][x][dim2 - 3][z] = -Txy[t1][x][dim2 - 4][z];
Txy[t1][x][dim2 - 2][z] = -Txy[t1][x][dim2 - 5][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyz[t1][x][1][z] = -Tyz[t1][x][2][z];
Tyz[t1][x][0][z] = -Tyz[t1][x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyz[t1][x][dim2 - 3][z] = -Tyz[t1][x][dim2 - 4][z];
Tyz[t1][x][dim2 - 2][z] = -Tyz[t1][x][dim2 - 5][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tyz[t1][x][y][1] = -Tyz[t1][x][y][2];
Tyz[t1][x][y][0] = -Tyz[t1][x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tyz[t1][x][y][dim3 - 3] = -Tyz[t1][x][y][dim3 - 4];
Tyz[t1][x][y][dim3 - 2] = -Tyz[t1][x][y][dim3 - 5];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txz[t1][1][y][z] = -Txz[t1][2][y][z];
Txz[t1][0][y][z] = -Txz[t1][3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txz[t1][dim1 - 3][y][z] = -Txz[t1][dim1 - 4][y][z];
Txz[t1][dim1 - 2][y][z] = -Txz[t1][dim1 - 5][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Txz[t1][x][y][1] = -Txz[t1][x][y][2];
Txz[t1][x][y][0] = -Txz[t1][x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Txz[t1][x][y][dim3 - 3] = -Txz[t1][x][y][dim3 - 4];
Txz[t1][x][y][dim3 - 2] = -Txz[t1][x][y][dim3 - 5];
}
}

#pragma omp for
for(int x=2;x<dim1 - 2;++x){
for(int y=2;y<dim2 - 2;++y){
#pragma ivdep
for(int z=2;z<dim3 - 2;++z){
U[t1][x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(27*Txz[t1][x][y][z] + Txz[t1][x][y][z - 2] - 27*Txz[t1][x][y][z - 1] - Txz[t1][x][y][z + 1]) + beta*dt*dx1*dx3*(27*Txy[t1][x][y][z] + Txy[t1][x][y - 2][z] - 27*Txy[t1][x][y - 1][z] - Txy[t1][x][y + 1][z]) + beta*dt*dx2*dx3*(-27*Txx[t1][x][y][z] + Txx[t1][x - 1][y][z] + 27*Txx[t1][x + 1][y][z] - Txx[t1][x + 2][y][z]) + 24*dx1*dx2*dx3*U[t0][x][y][z])/(dx1*dx2*dx3);
V[t1][x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(27*Tyz[t1][x][y][z] + Tyz[t1][x][y][z - 2] - 27*Tyz[t1][x][y][z - 1] - Tyz[t1][x][y][z + 1]) + beta*dt*dx1*dx3*(-27*Tyy[t1][x][y][z] + Tyy[t1][x][y - 1][z] + 27*Tyy[t1][x][y + 1][z] - Tyy[t1][x][y + 2][z]) + beta*dt*dx2*dx3*(27*Txy[t1][x][y][z] + Txy[t1][x - 2][y][z] - 27*Txy[t1][x - 1][y][z] - Txy[t1][x + 1][y][z]) + 24*dx1*dx2*dx3*V[t0][x][y][z])/(dx1*dx2*dx3);
W[t1][x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(-27*Tzz[t1][x][y][z] + Tzz[t1][x][y][z - 1] + 27*Tzz[t1][x][y][z + 1] - Tzz[t1][x][y][z + 2]) + beta*dt*dx1*dx3*(27*Tyz[t1][x][y][z] + Tyz[t1][x][y - 2][z] - 27*Tyz[t1][x][y - 1][z] - Tyz[t1][x][y + 1][z]) + beta*dt*dx2*dx3*(27*Txz[t1][x][y][z] + Txz[t1][x - 2][y][z] - 27*Txz[t1][x - 1][y][z] - Txz[t1][x + 1][y][z]) + 24*dx1*dx2*dx3*W[t0][x][y][z])/(dx1*dx2*dx3);
}
}
}

#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[t1][1][y][z] = (lambda*dx1*dx2*W[t1][2][y][z] - lambda*dx1*dx2*W[t1][2][y][z - 1] + lambda*dx1*dx3*V[t1][2][y][z] - lambda*dx1*dx3*V[t1][2][y - 1][z] + lambda*dx2*dx3*U[t1][2][y][z] + 2*mu*dx2*dx3*U[t1][2][y][z])/(dx2*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[t1][dim1 - 3][y][z] = (-lambda*dx1*dx2*W[t1][dim1 - 3][y][z] + lambda*dx1*dx2*W[t1][dim1 - 3][y][z - 1] - lambda*dx1*dx3*V[t1][dim1 - 3][y][z] + lambda*dx1*dx3*V[t1][dim1 - 3][y - 1][z] + lambda*dx2*dx3*U[t1][dim1 - 4][y][z] + 2*mu*dx2*dx3*U[t1][dim1 - 4][y][z])/(dx2*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[t1][1][y][z] = (-dx1*U[t1][1][y][z] + dx1*U[t1][1][y + 1][z] + dx1*U[t1][2][y][z] - dx1*U[t1][2][y + 1][z] + dx2*(2*V[t1][2][y][z] - V[t1][3][y][z]))/dx2;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[t1][dim1 - 2][y][z] = (-dx1*U[t1][dim1 - 4][y][z] + dx1*U[t1][dim1 - 4][y + 1][z] + dx1*U[t1][dim1 - 3][y][z] - dx1*U[t1][dim1 - 3][y + 1][z] + dx2*(-V[t1][dim1 - 4][y][z] + 2*V[t1][dim1 - 3][y][z]))/dx2;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[t1][1][y][z] = (-dx1*U[t1][1][y][z] + dx1*U[t1][1][y][z + 1] + dx1*U[t1][2][y][z] - dx1*U[t1][2][y][z + 1] + dx3*(2*W[t1][2][y][z] - W[t1][3][y][z]))/dx3;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[t1][dim1 - 2][y][z] = (-dx1*U[t1][dim1 - 4][y][z] + dx1*U[t1][dim1 - 4][y][z + 1] + dx1*U[t1][dim1 - 3][y][z] - dx1*U[t1][dim1 - 3][y][z + 1] + dx3*(-W[t1][dim1 - 4][y][z] + 2*W[t1][dim1 - 3][y][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[t1][x][1][z] = (lambda*dx1*dx2*W[t1][x][2][z] - lambda*dx1*dx2*W[t1][x][2][z - 1] + lambda*dx1*dx3*V[t1][x][2][z] + lambda*dx2*dx3*U[t1][x][2][z] - lambda*dx2*dx3*U[t1][x - 1][2][z] + 2*mu*dx1*dx3*V[t1][x][2][z])/(dx1*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[t1][x][dim2 - 3][z] = (-lambda*dx1*dx2*W[t1][x][dim2 - 3][z] + lambda*dx1*dx2*W[t1][x][dim2 - 3][z - 1] + lambda*dx1*dx3*V[t1][x][dim2 - 4][z] - lambda*dx2*dx3*U[t1][x][dim2 - 3][z] + lambda*dx2*dx3*U[t1][x - 1][dim2 - 3][z] + 2*mu*dx1*dx3*V[t1][x][dim2 - 4][z])/(dx1*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[t1][x][1][z] = (dx1*(2*U[t1][x][2][z] - U[t1][x][3][z]) - dx2*V[t1][x][1][z] + dx2*V[t1][x][2][z] + dx2*V[t1][x + 1][1][z] - dx2*V[t1][x + 1][2][z])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[t1][x][dim2 - 2][z] = (dx1*(-U[t1][x][dim2 - 4][z] + 2*U[t1][x][dim2 - 3][z]) - dx2*V[t1][x][dim2 - 4][z] + dx2*V[t1][x][dim2 - 3][z] + dx2*V[t1][x + 1][dim2 - 4][z] - dx2*V[t1][x + 1][dim2 - 3][z])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[t1][x][1][z] = (-dx2*V[t1][x][1][z] + dx2*V[t1][x][1][z + 1] + dx2*V[t1][x][2][z] - dx2*V[t1][x][2][z + 1] + dx3*(2*W[t1][x][2][z] - W[t1][x][3][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[t1][x][dim2 - 2][z] = (-dx2*V[t1][x][dim2 - 4][z] + dx2*V[t1][x][dim2 - 4][z + 1] + dx2*V[t1][x][dim2 - 3][z] - dx2*V[t1][x][dim2 - 3][z + 1] + dx3*(-W[t1][x][dim2 - 4][z] + 2*W[t1][x][dim2 - 3][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
W[t1][x][y][1] = (lambda*dx1*dx2*W[t1][x][y][2] + lambda*dx1*dx3*V[t1][x][y][2] - lambda*dx1*dx3*V[t1][x][y - 1][2] + lambda*dx2*dx3*U[t1][x][y][2] - lambda*dx2*dx3*U[t1][x - 1][y][2] + 2*mu*dx1*dx2*W[t1][x][y][2])/(dx1*dx2*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
W[t1][x][y][dim3 - 3] = (lambda*dx1*dx2*W[t1][x][y][dim3 - 4] - lambda*dx1*dx3*V[t1][x][y][dim3 - 3] + lambda*dx1*dx3*V[t1][x][y - 1][dim3 - 3] - lambda*dx2*dx3*U[t1][x][y][dim3 - 3] + lambda*dx2*dx3*U[t1][x - 1][y][dim3 - 3] + 2*mu*dx1*dx2*W[t1][x][y][dim3 - 4])/(dx1*dx2*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
U[t1][x][y][1] = (dx1*(2*U[t1][x][y][2] - U[t1][x][y][3]) - dx3*W[t1][x][y][1] + dx3*W[t1][x][y][2] + dx3*W[t1][x + 1][y][1] - dx3*W[t1][x + 1][y][2])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
U[t1][x][y][dim3 - 2] = (dx1*(-U[t1][x][y][dim3 - 4] + 2*U[t1][x][y][dim3 - 3]) - dx3*W[t1][x][y][dim3 - 4] + dx3*W[t1][x][y][dim3 - 3] + dx3*W[t1][x + 1][y][dim3 - 4] - dx3*W[t1][x + 1][y][dim3 - 3])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
V[t1][x][y][1] = (dx2*(2*V[t1][x][y][2] - V[t1][x][y][3]) - dx3*W[t1][x][y][1] + dx3*W[t1][x][y][2] + dx3*W[t1][x][y + 1][1] - dx3*W[t1][x][y + 1][2])/dx2;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
V[t1][x][y][dim3 - 2] = (dx2*(-V[t1][x][y][dim3 - 4] + 2*V[t1][x][y][dim3 - 3]) - dx3*W[t1][x][y][dim3 - 4] + dx3*W[t1][x][y][dim3 - 3] + dx3*W[t1][x][y + 1][dim3 - 4] - dx3*W[t1][x][y + 1][dim3 - 3])/dx2;
}
}




} // end of time loop
} // end of parallel section

printf("0.01\n");
printf("0.01\n");
printf("0.01\n");
float Txx_l2 = 0.0;
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Txx_l2+=pow(sqrt(2)*sqrt(mu/beta)*(sin(M_PI*y) - sin(M_PI*z))*sin(M_PI*x)*sin(sqrt(2)*M_PI*sqrt(beta*mu)) + Txx[0][_x][_y][_z], 2.0);
}
}
}
printf("Txx_l2\t%.10f\n", pow(Txx_l2*1.00000000000000e-6, 0.5));
float Tyy_l2 = 0.0;
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Tyy_l2+=pow(sqrt(2)*sqrt(mu/beta)*(-sin(M_PI*x) + sin(M_PI*z))*sin(M_PI*y)*sin(sqrt(2)*M_PI*sqrt(beta*mu)) + Tyy[0][_x][_y][_z], 2.0);
}
}
}
printf("Tyy_l2\t%.10f\n", pow(Tyy_l2*1.00000000000000e-6, 0.5));
float Tzz_l2 = 0.0;
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Tzz_l2+=pow(sqrt(2)*sqrt(mu/beta)*(sin(M_PI*x) - sin(M_PI*y))*sin(M_PI*z)*sin(sqrt(2)*M_PI*sqrt(beta*mu)) + Tzz[0][_x][_y][_z], 2.0);
}
}
}
printf("Tzz_l2\t%.10f\n", pow(Tzz_l2*1.00000000000000e-6, 0.5));
float Txy_l2 = 0.0;
for(int _x=2;_x<dim1 - 3;++_x){
float x= dx1*(_x - 1.5);
for(int _y=2;_y<dim2 - 3;++_y){
float y= dx2*(_y - 1.5);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Txy_l2+=pow(Txy[0][_x][_y][_z], 2.0);
}
}
}
printf("Txy_l2\t%.10f\n", pow(Txy_l2*1.00000000000000e-6, 0.5));
float Tyz_l2 = 0.0;
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 3;++_y){
float y= dx2*(_y - 1.5);
for(int _z=2;_z<dim3 - 3;++_z){
float z= dx3*(_z - 1.5);
Tyz_l2+=pow(Tyz[0][_x][_y][_z], 2.0);
}
}
}
printf("Tyz_l2\t%.10f\n", pow(Tyz_l2*1.00000000000000e-6, 0.5));
float Txz_l2 = 0.0;
for(int _x=2;_x<dim1 - 3;++_x){
float x= dx1*(_x - 1.5);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 3;++_z){
float z= dx3*(_z - 1.5);
Txz_l2+=pow(Txz[0][_x][_y][_z], 2.0);
}
}
}
printf("Txz_l2\t%.10f\n", pow(Txz_l2*1.00000000000000e-6, 0.5));
float U_l2 = 0.0;
for(int _x=2;_x<dim1 - 3;++_x){
float x= dx1*(_x - 1.5);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
U_l2+=pow(-(sin(M_PI*y) - sin(M_PI*z))*cos(M_PI*x)*cos(1.001*sqrt(2)*M_PI*sqrt(beta*mu)) + U[0][_x][_y][_z], 2.0);
}
}
}
printf("U_l2\t%.10f\n", pow(U_l2*1.00000000000000e-6, 0.5));
float V_l2 = 0.0;
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 3;++_y){
float y= dx2*(_y - 1.5);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
V_l2+=pow(-(-sin(M_PI*x) + sin(M_PI*z))*cos(M_PI*y)*cos(1.001*sqrt(2)*M_PI*sqrt(beta*mu)) + V[0][_x][_y][_z], 2.0);
}
}
}
printf("V_l2\t%.10f\n", pow(V_l2*1.00000000000000e-6, 0.5));
float W_l2 = 0.0;
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 3;++_z){
float z= dx3*(_z - 1.5);
W_l2+=pow(-(sin(M_PI*x) - sin(M_PI*y))*cos(M_PI*z)*cos(1.001*sqrt(2)*M_PI*sqrt(beta*mu)) + W[0][_x][_y][_z], 2.0);
}
}
}
printf("W_l2\t%.10f\n", pow(W_l2*1.00000000000000e-6, 0.5));


return 0;
}