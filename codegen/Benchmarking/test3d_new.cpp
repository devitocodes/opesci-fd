
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


//#define M_PI 3.14159265358979323846
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
const int tp = 2;
const float dt = 0.002;
const int margin = 2;
const int ntsteps = 1000;
const int vec_size = dim1*dim2*dim3;
const float mu = 0.5;
const float beta = 1.0;
const float rho = 1.0;
const float lambda = 0.5;
int t0 = 0;
int t1 = 0;

std::vector<float> _Txx_vec_1(vec_size);
float (*Txx_1)[dim2][dim3]= (float (*)[dim2][dim3]) _Txx_vec_1.data();
std::vector<float> _Tyy_vec_1(vec_size);
float (*Tyy_1)[dim2][dim3]= (float (*)[dim2][dim3]) _Tyy_vec_1.data();
std::vector<float> _Tzz_vec_1(vec_size);
float (*Tzz_1)[dim2][dim3]= (float (*)[dim2][dim3]) _Tzz_vec_1.data();
std::vector<float> _Txy_vec_1(vec_size);
float (*Txy_1)[dim2][dim3]= (float (*)[dim2][dim3]) _Txy_vec_1.data();
std::vector<float> _Tyz_vec_1(vec_size);
float (*Tyz_1)[dim2][dim3]= (float (*)[dim2][dim3]) _Tyz_vec_1.data();
std::vector<float> _Txz_vec_1(vec_size);
float (*Txz_1)[dim2][dim3]= (float (*)[dim2][dim3]) _Txz_vec_1.data();
std::vector<float> _U_vec_1(vec_size);
float (*U_1)[dim2][dim3]= (float (*)[dim2][dim3]) _U_vec_1.data();
std::vector<float> _V_vec_1(vec_size);
float (*V_1)[dim2][dim3]= (float (*)[dim2][dim3]) _V_vec_1.data();
std::vector<float> _W_vec_1(vec_size);
float (*W_1)[dim2][dim3]= (float (*)[dim2][dim3]) _W_vec_1.data();

std::vector<float> _Txx_vec_2(vec_size);
float (*Txx_2)[dim2][dim3]= (float (*)[dim2][dim3]) _Txx_vec_2.data();
std::vector<float> _Tyy_vec_2(vec_size);
float (*Tyy_2)[dim2][dim3]= (float (*)[dim2][dim3]) _Tyy_vec_2.data();
std::vector<float> _Tzz_vec_2(vec_size);
float (*Tzz_2)[dim2][dim3]= (float (*)[dim2][dim3]) _Tzz_vec_2.data();
std::vector<float> _Txy_vec_2(vec_size);
float (*Txy_2)[dim2][dim3]= (float (*)[dim2][dim3]) _Txy_vec_2.data();
std::vector<float> _Tyz_vec_2(vec_size);
float (*Tyz_2)[dim2][dim3]= (float (*)[dim2][dim3]) _Tyz_vec_2.data();
std::vector<float> _Txz_vec_2(vec_size);
float (*Txz_2)[dim2][dim3]= (float (*)[dim2][dim3]) _Txz_vec_2.data();
std::vector<float> _U_vec_2(vec_size);
float (*U_2)[dim2][dim3]= (float (*)[dim2][dim3]) _U_vec_2.data();
std::vector<float> _V_vec_2(vec_size);
float (*V_2)[dim2][dim3]= (float (*)[dim2][dim3]) _V_vec_2.data();
std::vector<float> _W_vec_2(vec_size);
float (*W_2)[dim2][dim3]= (float (*)[dim2][dim3]) _W_vec_2.data();

float (*Txx)[dim2][dim3];
float (*Tyy)[dim2][dim3];
float (*Tzz)[dim2][dim3];
float (*Txy)[dim2][dim3];
float (*Tyz)[dim2][dim3];
float (*Txz)[dim2][dim3];

float (*Txx_old)[dim2][dim3];
float (*Tyy_old)[dim2][dim3];
float (*Tzz_old)[dim2][dim3];
float (*Txy_old)[dim2][dim3];
float (*Tyz_old)[dim2][dim3];
float (*Txz_old)[dim2][dim3];
float (*U)[dim2][dim3];
float (*V)[dim2][dim3];
float (*W)[dim2][dim3];


float (*U_old)[dim2][dim3];
float (*V_old)[dim2][dim3];
float (*W_old)[dim2][dim3];


#pragma omp parallel
  {
	   	Txx = Txx_1;
	   	Tyy = Tyy_1;
	   	Tzz = Tzz_1;
	   	Txy = Txy_1;
	   	Tyz = Tyz_1;
	   	Txz = Txz_1;

	   	U = U_1;
	   	V = V_1;
	   	W = W_1;

	  
  #pragma omp for
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
Txx[_x][_y][_z]=0;
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
Tyy[_x][_y][_z]=0;
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
Tzz[_x][_y][_z]=0;
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
Txy[_x][_y][_z]=0.0;
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
Tyz[_x][_y][_z]=0.0;
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
Txz[_x][_y][_z]=0.0;
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
U[_x][_y][_z]=(sin(M_PI*y) - sin(M_PI*z))*cos(M_PI*x)*cos(0.001*sqrt(2)*M_PI*sqrt(mu/rho));
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
V[_x][_y][_z]=(-sin(M_PI*x) + sin(M_PI*z))*cos(M_PI*y)*cos(0.001*sqrt(2)*M_PI*sqrt(mu/rho));
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
W[_x][_y][_z]=(sin(M_PI*x) - sin(M_PI*y))*cos(M_PI*z)*cos(0.001*sqrt(2)*M_PI*sqrt(mu/rho));
}
}
}

  #pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txx[2][y][z] = 0;
Txx[1][y][z] = -Txx[3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txx[dim1 - 3][y][z] = 0;
Txx[dim1 - 2][y][z] = -Txx[dim1 - 4][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyy[x][2][z] = 0;
Tyy[x][1][z] = -Tyy[x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyy[x][dim2 - 3][z] = 0;
Tyy[x][dim2 - 2][z] = -Tyy[x][dim2 - 4][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tzz[x][y][2] = 0;
Tzz[x][y][1] = -Tzz[x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tzz[x][y][dim3 - 3] = 0;
Tzz[x][y][dim3 - 2] = -Tzz[x][y][dim3 - 4];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[1][y][z] = -Txy[2][y][z];
Txy[0][y][z] = -Txy[3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[dim1 - 3][y][z] = -Txy[dim1 - 4][y][z];
Txy[dim1 - 2][y][z] = -Txy[dim1 - 5][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[x][1][z] = -Txy[x][2][z];
Txy[0][x][z] = -Txy[x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[x][dim2 - 3][z] = -Txy[x][dim2 - 4][z];
Txy[x][dim2 - 2][z] = -Txy[x][dim2 - 5][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyz[x][1][z] = -Tyz[x][2][z];
Tyz[0][x][z] = -Tyz[x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyz[x][dim2 - 3][z] = -Tyz[x][dim2 - 4][z];
Tyz[x][dim2 - 2][z] = -Tyz[x][dim2 - 5][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tyz[x][y][1] = -Tyz[x][y][2];
Tyz[0][x][y] = -Tyz[x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tyz[x][y][dim3 - 3] = -Tyz[x][y][dim3 - 4];
Tyz[x][y][dim3 - 2] = -Tyz[x][y][dim3 - 5];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txz[1][y][z] = -Txz[2][y][z];
Txz[0][y][z] = -Txz[3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txz[dim1 - 3][y][z] = -Txz[dim1 - 4][y][z];
Txz[dim1 - 2][y][z] = -Txz[dim1 - 5][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Txz[x][y][1] = -Txz[x][y][2];
Txz[0][x][y] = -Txz[x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Txz[x][y][dim3 - 3] = -Txz[x][y][dim3 - 4];
Txz[x][y][dim3 - 2] = -Txz[x][y][dim3 - 5];
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[1][y][z] = (lambda*dx1*dx2*W[2][y][z] - lambda*dx1*dx2*W[2][y][z - 1] + lambda*dx1*dx3*V[2][y][z] - lambda*dx1*dx3*V[2][y - 1][z] + lambda*dx2*dx3*U[2][y][z] + 2*mu*dx2*dx3*U[2][y][z])/(dx2*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[dim1 - 3][y][z] = (-lambda*dx1*dx2*W[dim1 - 3][y][z] + lambda*dx1*dx2*W[dim1 - 3][y][z - 1] - lambda*dx1*dx3*V[dim1 - 3][y][z] + lambda*dx1*dx3*V[dim1 - 3][y - 1][z] + lambda*dx2*dx3*U[dim1 - 4][y][z] + 2*mu*dx2*dx3*U[dim1 - 4][y][z])/(dx2*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[1][y][z] = (-dx1*U[1][y][z] + dx1*U[1][y + 1][z] + dx1*U[2][y][z] - dx1*U[2][y + 1][z] + dx2*(2*V[2][y][z] - V[3][y][z]))/dx2;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[dim1 - 2][y][z] = (-dx1*U[dim1 - 4][y][z] + dx1*U[dim1 - 4][y + 1][z] + dx1*U[dim1 - 3][y][z] - dx1*U[dim1 - 3][y + 1][z] + dx2*(-V[dim1 - 4][y][z] + 2*V[dim1 - 3][y][z]))/dx2;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[1][y][z] = (-dx1*U[1][y][z] + dx1*U[1][y][z + 1] + dx1*U[2][y][z] - dx1*U[2][y][z + 1] + dx3*(2*W[2][y][z] - W[3][y][z]))/dx3;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[dim1 - 2][y][z] = (-dx1*U[dim1 - 4][y][z] + dx1*U[dim1 - 4][y][z + 1] + dx1*U[dim1 - 3][y][z] - dx1*U[dim1 - 3][y][z + 1] + dx3*(-W[dim1 - 4][y][z] + 2*W[dim1 - 3][y][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[x][1][z] = (lambda*dx1*dx2*W[x][2][z] - lambda*dx1*dx2*W[x][2][z - 1] + lambda*dx1*dx3*V[x][2][z] + lambda*dx2*dx3*U[x][2][z] - lambda*dx2*dx3*U[x - 1][2][z] + 2*mu*dx1*dx3*V[x][2][z])/(dx1*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[x][dim2 - 3][z] = (-lambda*dx1*dx2*W[x][dim2 - 3][z] + lambda*dx1*dx2*W[x][dim2 - 3][z - 1] + lambda*dx1*dx3*V[x][dim2 - 4][z] - lambda*dx2*dx3*U[x][dim2 - 3][z] + lambda*dx2*dx3*U[x - 1][dim2 - 3][z] + 2*mu*dx1*dx3*V[x][dim2 - 4][z])/(dx1*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[x][1][z] = (dx1*(2*U[x][2][z] - U[x][3][z]) - dx2*V[x][1][z] + dx2*V[x][2][z] + dx2*V[x + 1][1][z] - dx2*V[x + 1][2][z])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[x][dim2 - 2][z] = (dx1*(-U[x][dim2 - 4][z] + 2*U[x][dim2 - 3][z]) - dx2*V[x][dim2 - 4][z] + dx2*V[x][dim2 - 3][z] + dx2*V[x + 1][dim2 - 4][z] - dx2*V[x + 1][dim2 - 3][z])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[x][1][z] = (-dx2*V[x][1][z] + dx2*V[x][1][z + 1] + dx2*V[x][2][z] - dx2*V[x][2][z + 1] + dx3*(2*W[x][2][z] - W[x][3][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[x][dim2 - 2][z] = (-dx2*V[x][dim2 - 4][z] + dx2*V[x][dim2 - 4][z + 1] + dx2*V[x][dim2 - 3][z] - dx2*V[x][dim2 - 3][z + 1] + dx3*(-W[x][dim2 - 4][z] + 2*W[x][dim2 - 3][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
W[x][y][1] = (lambda*dx1*dx2*W[x][y][2] + lambda*dx1*dx3*V[x][y][2] - lambda*dx1*dx3*V[x][y - 1][2] + lambda*dx2*dx3*U[x][y][2] - lambda*dx2*dx3*U[x - 1][y][2] + 2*mu*dx1*dx2*W[x][y][2])/(dx1*dx2*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
W[x][y][dim3 - 3] = (lambda*dx1*dx2*W[x][y][dim3 - 4] - lambda*dx1*dx3*V[x][y][dim3 - 3] + lambda*dx1*dx3*V[x][y - 1][dim3 - 3] - lambda*dx2*dx3*U[x][y][dim3 - 3] + lambda*dx2*dx3*U[x - 1][y][dim3 - 3] + 2*mu*dx1*dx2*W[x][y][dim3 - 4])/(dx1*dx2*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
U[x][y][1] = (dx1*(2*U[x][y][2] - U[x][y][3]) - dx3*W[x][y][1] + dx3*W[x][y][2] + dx3*W[x + 1][y][1] - dx3*W[x + 1][y][2])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
U[x][y][dim3 - 2] = (dx1*(-U[x][y][dim3 - 4] + 2*U[x][y][dim3 - 3]) - dx3*W[x][y][dim3 - 4] + dx3*W[x][y][dim3 - 3] + dx3*W[x + 1][y][dim3 - 4] - dx3*W[x + 1][y][dim3 - 3])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
V[x][y][1] = (dx2*(2*V[x][y][2] - V[x][y][3]) - dx3*W[x][y][1] + dx3*W[x][y][2] + dx3*W[x][y + 1][1] - dx3*W[x][y + 1][2])/dx2;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
V[x][y][dim3 - 2] = (dx2*(-V[x][y][dim3 - 4] + 2*V[x][y][dim3 - 3]) - dx3*W[x][y][dim3 - 4] + dx3*W[x][y][dim3 - 3] + dx3*W[x][y + 1][dim3 - 4] - dx3*W[x][y + 1][dim3 - 3])/dx2;
}
}

// before this only touched txx_1 etc. 

  for(int _ti=0;_ti<ntsteps;_ti++){
	   if(_ti%2==0){
	   	Txx = Txx_2;
	   	Tyy = Tyy_2;
	   	Tzz = Tzz_2;
	   	Txy = Txy_2;
	   	Tyz = Tyz_2;
	   	Txz = Txz_2;

	   	Txx_old = Txx_1;
	   	Tyy_old = Tyy_1;
	   	Tzz_old = Tzz_1;
	   	Txy_old = Txy_1;
	   	Tyz_old = Tyz_1;
	   	Txz_old = Txz_1;

	   	U = U_1;
	   	V = V_1;
	   	W = W_1;

		U_old = U_2;
	   	V_old = V_2;
	   	W_old = W_2;

	   }else{
		Txx = Txx_1;
	   	Tyy = Tyy_1;
	   	Tzz = Tzz_1;
	   	Txy = Txy_1;
	   	Tyz = Tyz_1;
	   	Txz = Txz_1;

	   	Txx_old = Txx_2;
	   	Tyy_old = Tyy_2;
	   	Tzz_old = Tzz_2;
	   	Txy_old = Txy_2;
	   	Tyz_old = Tyz_2;
	   	Txz_old = Txz_2;

	   	U = U_2;
	   	V = V_2;
	   	W = W_2;

	   	U_old = U_1;
	   	V_old = V_1;
	   	W_old = W_1;

	   }

//critical(0,t0,t1);

    #pragma omp for
for(int x=2;x<dim1 - 2;++x){
for(int y=2;y<dim2 - 2;++y){
#pragma ivdep
for(int z=2;z<dim3 - 2;++z){
//Txx[x][y][z] = U[x][y][z+1];
Txx[x][y][z]= (1.0F/24.0F)*(lambda*dt*dx1*dx2*(27*W[x][y][z] + W[x][y][z - 2] - 27*W[x][y][z - 1] - W[x][y][z + 1]) + lambda*dt*dx1*dx3*(27*V[x][y][z] + V[x][y - 2][z] - 27*V[x][y - 1][z] - V[x][y + 1][z]) + dt*dx2*dx3*(27*lambda*U[x][y][z] + lambda*U[x - 2][y][z] - 27*lambda*U[x - 1][y][z] - lambda*U[x + 1][y][z] + 54*mu*U[x][y][z] + 2*mu*U[x - 2][y][z] - 54*mu*U[x - 1][y][z] - 2*mu*U[x + 1][y][z]) + 24*dx1*dx2*dx3*Txx_old[x][y][z])/(dx1*dx2*dx3);
Tyy[x][y][z]=(1.0F/24.0F)*(lambda*dt*dx1*dx2*(27*W[x][y][z] + W[x][y][z - 2] - 27*W[x][y][z - 1] - W[x][y][z + 1]) + lambda*dt*dx2*dx3*(27*U[x][y][z] + U[x - 2][y][z] - 27*U[x - 1][y][z] - U[x + 1][y][z]) + dt*dx1*dx3*(27*lambda*V[x][y][z] + lambda*V[x][y - 2][z] - 27*lambda*V[x][y - 1][z] - lambda*V[x][y + 1][z] + 54*mu*V[x][y][z] + 2*mu*V[x][y - 2][z] - 54*mu*V[x][y - 1][z] - 2*mu*V[x][y + 1][z]) + 24*dx1*dx2*dx3*Tyy_old[x][y][z])/(dx1*dx2*dx3);
Tzz[x][y][z]=(1.0F/24.0F)*(lambda*dt*dx1*dx3*(27*V[x][y][z] + V[x][y - 2][z] - 27*V[x][y - 1][z] - V[x][y + 1][z]) + lambda*dt*dx2*dx3*(27*U[x][y][z] + U[x - 2][y][z] - 27*U[x - 1][y][z] - U[x + 1][y][z]) + dt*dx1*dx2*(27*lambda*W[x][y][z] + lambda*W[x][y][z - 2] - 27*lambda*W[x][y][z - 1] - lambda*W[x][y][z + 1] + 54*mu*W[x][y][z] + 2*mu*W[x][y][z - 2] - 54*mu*W[x][y][z - 1] - 2*mu*W[x][y][z + 1]) + 24*dx1*dx2*dx3*Tzz_old[x][y][z])/(dx1*dx2*dx3);
Txy[x][y][z]=(1.0F/24.0F)*(mu*dt*dx1*(-27*U[x][y][z] + U[x][y - 1][z] + 27*U[x][y + 1][z] - U[x][y + 2][z]) + mu*dt*dx2*(-27*V[x][y][z] + V[x - 1][y][z] + 27*V[x + 1][y][z] - V[x + 2][y][z]) + 24*dx1*dx2*Txy_old[x][y][z])/(dx1*dx2);
Tyz[x][y][z]=(1.0F/24.0F)*(mu*dt*dx2*(-27*V[x][y][z] + V[x][y][z - 1] + 27*V[x][y][z + 1] - V[x][y][z + 2]) + mu*dt*dx3*(-27*W[x][y][z] + W[x][y - 1][z] + 27*W[x][y + 1][z] - W[x][y + 2][z]) + 24*dx2*dx3*Tyz_old[x][y][z])/(dx2*dx3);
Txz[x][y][z]=(1.0F/24.0F)*(mu*dt*dx1*(-27*U[x][y][z] + U[x][y][z - 1] + 27*U[x][y][z + 1] - U[x][y][z + 2]) + mu*dt*dx3*(-27*W[x][y][z] + W[x - 1][y][z] + 27*W[x + 1][y][z] - W[x + 2][y][z]) + 24*dx1*dx3*Txz_old[x][y][z])/(dx1*dx3);
}
}
}

/*#pragma omp for
for (int c0 = 0; c0 <= 3; c0 += 1){
  for (int c1 = 0; c1 <= 3; c1 += 1){
    for (int c2 = 0; c2 <= 3; c2 += 1){
      for (int c3 = 0; c3 <= fmin(31, -32 * c0 + 100); c3 += 1){
        for (int c4 = 0; c4 <= fmin(31, -32 * c1 + 100); c4 += 1){
          for (int c5 = 0; c5 <= fmin(31, -32 * c2 + 100); c5 += 1){
          	int x = 32 * c0 + c3 + 2 ;
          	int y =  32 * c1 + c4 + 2;
          	int z =  32 * c2 + c5 + 2;
          	Txx[x][y][z]= (1.0F/24.0F)*(lambda*dt*dx1*dx2*(27*W[x][y][z] + W[x][y][z - 2] - 27*W[x][y][z - 1] - W[x][y][z + 1]) + lambda*dt*dx1*dx3*(27*V[x][y][z] + V[x][y - 2][z] - 27*V[x][y - 1][z] - V[x][y + 1][z]) + dt*dx2*dx3*(27*lambda*U[x][y][z] + lambda*U[x - 2][y][z] - 27*lambda*U[x - 1][y][z] - lambda*U[x + 1][y][z] + 54*mu*U[x][y][z] + 2*mu*U[x - 2][y][z] - 54*mu*U[x - 1][y][z] - 2*mu*U[x + 1][y][z]) + 24*dx1*dx2*dx3*Txx_old[x][y][z])/(dx1*dx2*dx3);
			Tyy[x][y][z]=(1.0F/24.0F)*(lambda*dt*dx1*dx2*(27*W[x][y][z] + W[x][y][z - 2] - 27*W[x][y][z - 1] - W[x][y][z + 1]) + lambda*dt*dx2*dx3*(27*U[x][y][z] + U[x - 2][y][z] - 27*U[x - 1][y][z] - U[x + 1][y][z]) + dt*dx1*dx3*(27*lambda*V[x][y][z] + lambda*V[x][y - 2][z] - 27*lambda*V[x][y - 1][z] - lambda*V[x][y + 1][z] + 54*mu*V[x][y][z] + 2*mu*V[x][y - 2][z] - 54*mu*V[x][y - 1][z] - 2*mu*V[x][y + 1][z]) + 24*dx1*dx2*dx3*Tyy_old[x][y][z])/(dx1*dx2*dx3);
			Tzz[x][y][z]=(1.0F/24.0F)*(lambda*dt*dx1*dx3*(27*V[x][y][z] + V[x][y - 2][z] - 27*V[x][y - 1][z] - V[x][y + 1][z]) + lambda*dt*dx2*dx3*(27*U[x][y][z] + U[x - 2][y][z] - 27*U[x - 1][y][z] - U[x + 1][y][z]) + dt*dx1*dx2*(27*lambda*W[x][y][z] + lambda*W[x][y][z - 2] - 27*lambda*W[x][y][z - 1] - lambda*W[x][y][z + 1] + 54*mu*W[x][y][z] + 2*mu*W[x][y][z - 2] - 54*mu*W[x][y][z - 1] - 2*mu*W[x][y][z + 1]) + 24*dx1*dx2*dx3*Tzz_old[x][y][z])/(dx1*dx2*dx3);
			Txy[x][y][z]=(1.0F/24.0F)*(mu*dt*dx1*(-27*U[x][y][z] + U[x][y - 1][z] + 27*U[x][y + 1][z] - U[x][y + 2][z]) + mu*dt*dx2*(-27*V[x][y][z] + V[x - 1][y][z] + 27*V[x + 1][y][z] - V[x + 2][y][z]) + 24*dx1*dx2*Txy_old[x][y][z])/(dx1*dx2);
			Tyz[x][y][z]=(1.0F/24.0F)*(mu*dt*dx2*(-27*V[x][y][z] + V[x][y][z - 1] + 27*V[x][y][z + 1] - V[x][y][z + 2]) + mu*dt*dx3*(-27*W[x][y][z] + W[x][y - 1][z] + 27*W[x][y + 1][z] - W[x][y + 2][z]) + 24*dx2*dx3*Tyz_old[x][y][z])/(dx2*dx3);
			Txz[x][y][z]=(1.0F/24.0F)*(mu*dt*dx1*(-27*U[x][y][z] + U[x][y][z - 1] + 27*U[x][y][z + 1] - U[x][y][z + 2]) + mu*dt*dx3*(-27*W[x][y][z] + W[x - 1][y][z] + 27*W[x + 1][y][z] - W[x + 2][y][z]) + 24*dx1*dx3*Txz_old[x][y][z])/(dx1*dx3);

           }
        }
      }
    }
  }
}
*/
if(_ti%2==0){
Txx = Txx_2;
Tyy = Tyy_2;
Tzz = Tzz_2;
Txy = Txy_2;
Tyz = Tyz_2;
Txz = Txz_2;
}else{
Txx = Txx_1;
Tyy = Tyy_1;
Tzz = Tzz_1;
Txy = Txy_1;
Tyz = Tyz_1;
Txz = Txz_1;

}


// before next computation only touched Txx ~ Txz [t1]
    #pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txx[2][y][z] = 0;
Txx[1][y][z] = -Txx[3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txx[dim1 - 3][y][z] = 0;
Txx[dim1 - 2][y][z] = -Txx[dim1 - 4][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyy[x][2][z] = 0;
Tyy[x][1][z] = -Tyy[x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyy[x][dim2 - 3][z] = 0;
Tyy[x][dim2 - 2][z] = -Tyy[x][dim2 - 4][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tzz[x][y][2] = 0;
Tzz[x][y][1] = -Tzz[x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tzz[x][y][dim3 - 3] = 0;
Tzz[x][y][dim3 - 2] = -Tzz[x][y][dim3 - 4];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[1][y][z] = -Txy[2][y][z];
Txy[0][y][z] = -Txy[3][y][z];
}
}

// TODO: below breaks tiling
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[dim1 - 3][y][z] = -Txy[dim1 - 4][y][z];
Txy[dim1 - 2][y][z] = -Txy[dim1 - 5][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[x][1][z] = -Txy[x][2][z];
Txy[x][0][z] = -Txy[x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txy[x][dim2 - 3][z] = -Txy[x][dim2 - 4][z];
Txy[x][dim2 - 2][z] = -Txy[x][dim2 - 5][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyz[x][1][z] = -Tyz[x][2][z];
Tyz[x][0][z] = -Tyz[x][3][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int z=0;z<dim3;++z){
Tyz[x][dim2 - 3][z] = -Tyz[x][dim2 - 4][z];
Tyz[x][dim2 - 2][z] = -Tyz[x][dim2 - 5][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tyz[x][y][1] = -Tyz[x][y][2];
Tyz[x][y][0] = -Tyz[x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Tyz[x][y][dim3 - 3] = -Tyz[x][y][dim3 - 4];
Tyz[x][y][dim3 - 2] = -Tyz[x][y][dim3 - 5];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txz[1][y][z] = -Txz[2][y][z];
Txz[0][y][z] = -Txz[3][y][z];
}
}
#pragma omp for
for(int y=0;y<dim2;++y){
#pragma ivdep
for(int z=0;z<dim3;++z){
Txz[dim1 - 3][y][z] = -Txz[dim1 - 4][y][z];
Txz[dim1 - 2][y][z] = -Txz[dim1 - 5][y][z];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Txz[x][y][1] = -Txz[x][y][2];
Txz[x][y][0] = -Txz[x][y][3];
}
}
#pragma omp for
for(int x=0;x<dim1;++x){
#pragma ivdep
for(int y=0;y<dim2;++y){
Txz[x][y][dim3 - 3] = -Txz[x][y][dim3 - 4];
Txz[x][y][dim3 - 2] = -Txz[x][y][dim3 - 5];
}
}
//critical(1,t0,t1);
//above works

if(_ti%2==0){
	   	Txx = Txx_2;
	   	Tyy = Tyy_2;
	   	Tzz = Tzz_2;
	   	Txy = Txy_2;
	   	Tyz = Tyz_2;
	   	Txz = Txz_2;
	   	U = U_2;
	   	V = V_2;
	   	W = W_2;

		U_old = U_1;
	   	V_old = V_1;
	   	W_old = W_1;

	   }else{
		Txx = Txx_1;
	   	Tyy = Tyy_1;
	   	Tzz = Tzz_1;
	   	Txy = Txy_1;
	   	Tyz = Tyz_1;
	   	Txz = Txz_1;
	   	U = U_1;
	   	V = V_1;
	   	W = W_1;

	   	U_old = U_2;
	   	V_old = V_2;
	   	W_old = W_2;

	   }

    #pragma omp for
for(int x=2;x<dim1 - 2;++x){
for(int y=2;y<dim2 - 2;++y){
#pragma ivdep
for(int z=2;z<dim3 - 2;++z){
U[x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(27*Txz[x][y][z] + Txz[x][y][z - 2] - 27*Txz[x][y][z - 1] - Txz[x][y][z + 1]) + beta*dt*dx1*dx3*(27*Txy[x][y][z] + Txy[x][y - 2][z] - 27*Txy[x][y - 1][z] - Txy[x][y + 1][z]) + beta*dt*dx2*dx3*(-27*Txx[x][y][z] + Txx[x - 1][y][z] + 27*Txx[x + 1][y][z] - Txx[x + 2][y][z]) + 24*dx1*dx2*dx3*U_old[x][y][z])/(dx1*dx2*dx3);
V[x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(27*Tyz[x][y][z] + Tyz[x][y][z - 2] - 27*Tyz[x][y][z - 1] - Tyz[x][y][z + 1]) + beta*dt*dx1*dx3*(-27*Tyy[x][y][z] + Tyy[x][y - 1][z] + 27*Tyy[x][y + 1][z] - Tyy[x][y + 2][z]) + beta*dt*dx2*dx3*(27*Txy[x][y][z] + Txy[x - 2][y][z] - 27*Txy[x - 1][y][z] - Txy[x + 1][y][z]) + 24*dx1*dx2*dx3*V_old[x][y][z])/(dx1*dx2*dx3);
W[x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(-27*Tzz[x][y][z] + Tzz[x][y][z - 1] + 27*Tzz[x][y][z + 1] - Tzz[x][y][z + 2]) + beta*dt*dx1*dx3*(27*Tyz[x][y][z] + Tyz[x][y - 2][z] - 27*Tyz[x][y - 1][z] - Tyz[x][y + 1][z]) + beta*dt*dx2*dx3*(27*Txz[x][y][z] + Txz[x - 2][y][z] - 27*Txz[x - 1][y][z] - Txz[x + 1][y][z]) + 24*dx1*dx2*dx3*W_old[x][y][z])/(dx1*dx2*dx3);
}
}
}
/*
#pragma omp for
for (int c0 = 0; c0 <= 3; c0 += 1){
  for (int c1 = 0; c1 <= 3; c1 += 1){
    for (int c2 = 0; c2 <= 3; c2 += 1){
      for (int c3 = 0; c3 <= fmin(31, -32 * c0 + 100); c3 += 1){
        for (int c4 = 0; c4 <= fmin(31, -32 * c1 + 100); c4 += 1){
          for (int c5 = 0; c5 <= fmin(31, -32 * c2 + 100); c5 += 1){
          	int x = 32 * c0 + c3 + 2 ;
          	int y =  32 * c1 + c4 + 2;
          	int z =  32 * c2 + c5 + 2;
          	U[x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(27*Txz[x][y][z] + Txz[x][y][z - 2] - 27*Txz[x][y][z - 1] - Txz[x][y][z + 1]) + beta*dt*dx1*dx3*(27*Txy[x][y][z] + Txy[x][y - 2][z] - 27*Txy[x][y - 1][z] - Txy[x][y + 1][z]) + beta*dt*dx2*dx3*(-27*Txx[x][y][z] + Txx[x - 1][y][z] + 27*Txx[x + 1][y][z] - Txx[x + 2][y][z]) + 24*dx1*dx2*dx3*U_old[x][y][z])/(dx1*dx2*dx3);
			V[x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(27*Tyz[x][y][z] + Tyz[x][y][z - 2] - 27*Tyz[x][y][z - 1] - Tyz[x][y][z + 1]) + beta*dt*dx1*dx3*(-27*Tyy[x][y][z] + Tyy[x][y - 1][z] + 27*Tyy[x][y + 1][z] - Tyy[x][y + 2][z]) + beta*dt*dx2*dx3*(27*Txy[x][y][z] + Txy[x - 2][y][z] - 27*Txy[x - 1][y][z] - Txy[x + 1][y][z]) + 24*dx1*dx2*dx3*V_old[x][y][z])/(dx1*dx2*dx3);
			W[x][y][z]=(1.0F/24.0F)*(beta*dt*dx1*dx2*(-27*Tzz[x][y][z] + Tzz[x][y][z - 1] + 27*Tzz[x][y][z + 1] - Tzz[x][y][z + 2]) + beta*dt*dx1*dx3*(27*Tyz[x][y][z] + Tyz[x][y - 2][z] - 27*Tyz[x][y - 1][z] - Tyz[x][y + 1][z]) + beta*dt*dx2*dx3*(27*Txz[x][y][z] + Txz[x - 2][y][z] - 27*Txz[x - 1][y][z] - Txz[x + 1][y][z]) + 24*dx1*dx2*dx3*W_old[x][y][z])/(dx1*dx2*dx3);



           }
        }
      }
    }
  }
}
*/

if(_ti%2==0){
U = U_2;
V = V_2;
W = W_2;
}else{

U = U_1;
V = V_1;
W = W_1;

}

    #pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[1][y][z] = (lambda*dx1*dx2*W[2][y][z] - lambda*dx1*dx2*W[2][y][z - 1] + lambda*dx1*dx3*V[2][y][z] - lambda*dx1*dx3*V[2][y - 1][z] + lambda*dx2*dx3*U[2][y][z] + 2*mu*dx2*dx3*U[2][y][z])/(dx2*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[dim1 - 3][y][z] = (-lambda*dx1*dx2*W[dim1 - 3][y][z] + lambda*dx1*dx2*W[dim1 - 3][y][z - 1] - lambda*dx1*dx3*V[dim1 - 3][y][z] + lambda*dx1*dx3*V[dim1 - 3][y - 1][z] + lambda*dx2*dx3*U[dim1 - 4][y][z] + 2*mu*dx2*dx3*U[dim1 - 4][y][z])/(dx2*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[1][y][z] = (-dx1*U[1][y][z] + dx1*U[1][y + 1][z] + dx1*U[2][y][z] - dx1*U[2][y + 1][z] + dx2*(2*V[2][y][z] - V[3][y][z]))/dx2;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[dim1 - 2][y][z] = (-dx1*U[dim1 - 4][y][z] + dx1*U[dim1 - 4][y + 1][z] + dx1*U[dim1 - 3][y][z] - dx1*U[dim1 - 3][y + 1][z] + dx2*(-V[dim1 - 4][y][z] + 2*V[dim1 - 3][y][z]))/dx2;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[1][y][z] = (-dx1*U[1][y][z] + dx1*U[1][y][z + 1] + dx1*U[2][y][z] - dx1*U[2][y][z + 1] + dx3*(2*W[2][y][z] - W[3][y][z]))/dx3;
}
}
#pragma omp for
for(int y=1;y<dim2 - 1;++y){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[dim1 - 2][y][z] = (-dx1*U[dim1 - 4][y][z] + dx1*U[dim1 - 4][y][z + 1] + dx1*U[dim1 - 3][y][z] - dx1*U[dim1 - 3][y][z + 1] + dx3*(-W[dim1 - 4][y][z] + 2*W[dim1 - 3][y][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[x][1][z] = (lambda*dx1*dx2*W[x][2][z] - lambda*dx1*dx2*W[x][2][z - 1] + lambda*dx1*dx3*V[x][2][z] + lambda*dx2*dx3*U[x][2][z] - lambda*dx2*dx3*U[x - 1][2][z] + 2*mu*dx1*dx3*V[x][2][z])/(dx1*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
V[x][dim2 - 3][z] = (-lambda*dx1*dx2*W[x][dim2 - 3][z] + lambda*dx1*dx2*W[x][dim2 - 3][z - 1] + lambda*dx1*dx3*V[x][dim2 - 4][z] - lambda*dx2*dx3*U[x][dim2 - 3][z] + lambda*dx2*dx3*U[x - 1][dim2 - 3][z] + 2*mu*dx1*dx3*V[x][dim2 - 4][z])/(dx1*dx3*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[x][1][z] = (dx1*(2*U[x][2][z] - U[x][3][z]) - dx2*V[x][1][z] + dx2*V[x][2][z] + dx2*V[x + 1][1][z] - dx2*V[x + 1][2][z])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
U[x][dim2 - 2][z] = (dx1*(-U[x][dim2 - 4][z] + 2*U[x][dim2 - 3][z]) - dx2*V[x][dim2 - 4][z] + dx2*V[x][dim2 - 3][z] + dx2*V[x + 1][dim2 - 4][z] - dx2*V[x + 1][dim2 - 3][z])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[x][1][z] = (-dx2*V[x][1][z] + dx2*V[x][1][z + 1] + dx2*V[x][2][z] - dx2*V[x][2][z + 1] + dx3*(2*W[x][2][z] - W[x][3][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int z=1;z<dim3 - 1;++z){
W[x][dim2 - 2][z] = (-dx2*V[x][dim2 - 4][z] + dx2*V[x][dim2 - 4][z + 1] + dx2*V[x][dim2 - 3][z] - dx2*V[x][dim2 - 3][z + 1] + dx3*(-W[x][dim2 - 4][z] + 2*W[x][dim2 - 3][z]))/dx3;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
W[x][y][1] = (lambda*dx1*dx2*W[x][y][2] + lambda*dx1*dx3*V[x][y][2] - lambda*dx1*dx3*V[x][y - 1][2] + lambda*dx2*dx3*U[x][y][2] - lambda*dx2*dx3*U[x - 1][y][2] + 2*mu*dx1*dx2*W[x][y][2])/(dx1*dx2*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
W[x][y][dim3 - 3] = (lambda*dx1*dx2*W[x][y][dim3 - 4] - lambda*dx1*dx3*V[x][y][dim3 - 3] + lambda*dx1*dx3*V[x][y - 1][dim3 - 3] - lambda*dx2*dx3*U[x][y][dim3 - 3] + lambda*dx2*dx3*U[x - 1][y][dim3 - 3] + 2*mu*dx1*dx2*W[x][y][dim3 - 4])/(dx1*dx2*(lambda + 2*mu));
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
U[x][y][1] = (dx1*(2*U[x][y][2] - U[x][y][3]) - dx3*W[x][y][1] + dx3*W[x][y][2] + dx3*W[x + 1][y][1] - dx3*W[x + 1][y][2])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
U[x][y][dim3 - 2] = (dx1*(-U[x][y][dim3 - 4] + 2*U[x][y][dim3 - 3]) - dx3*W[x][y][dim3 - 4] + dx3*W[x][y][dim3 - 3] + dx3*W[x + 1][y][dim3 - 4] - dx3*W[x + 1][y][dim3 - 3])/dx1;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
V[x][y][1] = (dx2*(2*V[x][y][2] - V[x][y][3]) - dx3*W[x][y][1] + dx3*W[x][y][2] + dx3*W[x][y + 1][1] - dx3*W[x][y + 1][2])/dx2;
}
}
#pragma omp for
for(int x=1;x<dim1 - 1;++x){
#pragma ivdep
for(int y=1;y<dim2 - 1;++y){
V[x][y][dim3 - 2] = (dx2*(-V[x][y][dim3 - 4] + 2*V[x][y][dim3 - 3]) - dx3*W[x][y][dim3 - 4] + dx3*W[x][y][dim3 - 3] + dx3*W[x][y + 1][dim3 - 4] - dx3*W[x][y + 1][dim3 - 3])/dx2;
}
}


    
//printf("%f\n", Txx[2][2][2]);

  } // end of time loop
  } // end of parallel section

  /*printf("0.01\n");
printf("0.01\n");
printf("0.01\n");*/

Txx = Txx_1;
Tyy = Tyy_1;
Tzz = Tzz_1;
Txy = Txy_1;
Tyz = Tyz_1;
Txz = Txz_1;

U = U_1;
V = V_1;
W = W_1;

float Txx_l2 = 0.0;
for(int _x=2;_x<dim1 - 2;++_x){
float x= dx1*(_x - 2);
for(int _y=2;_y<dim2 - 2;++_y){
float y= dx2*(_y - 2);
for(int _z=2;_z<dim3 - 2;++_z){
float z= dx3*(_z - 2);
//printf("%f \n",Txx[_x][_y][_z] );
Txx_l2+=pow(sqrt(2)*sqrt(mu*rho)*(sin(M_PI*y) - sin(M_PI*z))*sin(M_PI*x)*sin(2.0*sqrt(2)*M_PI*sqrt(mu/rho)) + Txx[_x][_y][_z], 2.0);
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
Tyy_l2+=pow(sqrt(2)*sqrt(mu*rho)*(-sin(M_PI*x) + sin(M_PI*z))*sin(M_PI*y)*sin(2.0*sqrt(2)*M_PI*sqrt(mu/rho)) + Tyy[_x][_y][_z], 2.0);
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
Tzz_l2+=pow(sqrt(2)*sqrt(mu*rho)*(sin(M_PI*x) - sin(M_PI*y))*sin(M_PI*z)*sin(2.0*sqrt(2)*M_PI*sqrt(mu/rho)) + Tzz[_x][_y][_z], 2.0);
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
Txy_l2+=pow(Txy[_x][_y][_z], 2.0);
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
Tyz_l2+=pow(Tyz[_x][_y][_z], 2.0);
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
Txz_l2+=pow(Txz[_x][_y][_z], 2.0);
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
U_l2+=pow(-(sin(M_PI*y) - sin(M_PI*z))*cos(M_PI*x)*cos(2.001*sqrt(2)*M_PI*sqrt(mu/rho)) + U[_x][_y][_z], 2.0);
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
V_l2+=pow(-(-sin(M_PI*x) + sin(M_PI*z))*cos(M_PI*y)*cos(2.001*sqrt(2)*M_PI*sqrt(mu/rho)) + V[_x][_y][_z], 2.0);
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
W_l2+=pow(-(sin(M_PI*x) - sin(M_PI*y))*cos(M_PI*z)*cos(2.001*sqrt(2)*M_PI*sqrt(mu/rho)) + W[_x][_y][_z], 2.0);
}
}
}
printf("W_l2\t%.10f\n", pow(W_l2*1.00000000000000e-6, 0.5));

  return 0;
}