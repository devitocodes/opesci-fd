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

#include <pochoir.hpp>

Pochoir_Boundary_3D(fd_bv_3D, arr, t, i, j, k)
    return 0;
Pochoir_Boundary_End

int main(){
  std::string vpfile("VPhomogx200");       // Vp file in binary (in m/s)
  std::string vsfile("VShomogx200");       // Vs file: Vs file in binary (in m/s)
  std::string rhofile("RHOhomogx200");     // rhofile: densities file in binary (in kg/m**3)
  std::string geomsrc("geomx200.src");     // Geometry file with sources locations (in m)
  std::string geomrec("geomx200.rec");     // Geometry file with receivers locations (in m)
  int dimx=200, dimy=200, dimz=200;        // Model dimensions in number of nodes
  float h=25;                              // Node spacing in meters (hx=hy=hz=h)
  std::string xsrcname("xrick5.sep");      // Source files
  std::string ysrcname("yrick5.sep");
  std::string zsrcname("zrick5.sep"); 
  float sdt=0.004;                         // Source sample rate (in s).
  float maxt=1.0;                          // Maximum time to compute;
  std::string useis("seisu");              // Seismogram files for the x, y and z components of velocity.
  std::string vseis("seisv");
  std::string wseis("seisw");
  std::string pseis("seisp");              // Seismogram files for pressure.

  int index = 0;
  
  // Read Vp
  std::vector<float> vp;
  if(opesci_read_simple_binary(vpfile.c_str(), vp)){
    opesci_abort("Missing input file\n");
  }
  assert(vp.size()==dimx*dimy*dimz);
 
  // Read Vs
  std::vector<float> vs;
  if(opesci_read_simple_binary(vsfile.c_str(), vs)){
    opesci_abort("Missing input file\n");
  }
  assert(vs.size()==dimx*dimy*dimz);

  // Read rho (density).
  std::vector<float> rho;
  if(opesci_read_simple_binary(rhofile.c_str(), rho)){
    opesci_abort("Missing input file\n");
  }
  assert(rho.size()==dimx*dimy*dimz);
  std::vector<float> buoyancy(rho.size());

  // Calculate Lame constents.
  std::vector<float> mu, lam;
  opesci_calculate_lame_costants(vp, vs, rho, mu, lam);

  // Get sources.
  std::vector<float> coorsrc, xsrc, ysrc, zsrc;
  if(opesci_read_souces(geomsrc.c_str(), xsrcname.c_str(), ysrcname.c_str(), zsrcname.c_str(), coorsrc, xsrc, ysrc, zsrc)){
    opesci_abort("Cannot read sources.\n");
  }
  int nsrc = coorsrc.size()/3;

  // Get receivers.
  std::vector<float> coorrec;
  if(opesci_read_receivers(geomrec.c_str(), coorrec)){
    opesci_abort("Cannot read receivers.\n");
  }
  int nrec = coorrec.size()/3;

  float dt = opesci_calculate_dt(vp, h);
  int ntsteps = (int)(maxt/dt);
  
  // Resample sources if required.
  assert(nsrc==1);
  std::vector<float> resampled_src;
  opesci_resample_timeseries(xsrc, dt, sdt, resampled_src);
  xsrc.swap(resampled_src);

  opesci_resample_timeseries(ysrc, dt, sdt, resampled_src);
  ysrc.swap(resampled_src);

  opesci_resample_timeseries(zsrc, dt, sdt, resampled_src);
  zsrc.swap(resampled_src);

  int snt = xsrc.size()/nsrc;
  sdt = dt;

  float c0=9./8.;
  float c1=1./24.;
  
  Pochoir_Shape_3D fd_shape_3D[] = {
    {1,0,0,0}, 
    {0,1,0,0},
    {0,0,1,0},
    {0,0,0,1},
    {0,2,0,0},
    {0,0,2,0},
    {0,0,0,2},
    {0,-1,0,0},
    {0,0,-1,0},
    {0,0,0,-1},
    {0,-2,0,0},
    {0,0,-2,0},
    {0,0,0,-2}};
  
  // Prognostic fields.
  Pochoir_Array<float, 3> u(dimz, dimy, dimx), v(dimz, dimy, dimx), w(dimz, dimy, dimx),
    txx(dimz, dimy, dimx), tyy(dimz, dimy, dimx), tzz(dimz, dimy, dimx),
    tyz(dimz, dimy, dimx), txz(dimz, dimy, dimx), txy(dimz, dimy, dimx);

  // Subsurface model
  Pochoir_Array<float, 3> Buoyancy(dimz, dimy, dimx), Lambda(dimz, dimy, dimx), Mu(dimz, dimy, dimx);
  
  Pochoir_3D fd_3D(fd_shape_3D);
  int ds=2;
  Pochoir_Domain I(0+ds, dimx-ds), J(0+ds, dimy-ds), K(0+ds, dimz-ds);
  
  u.Register_Boundary(fd_bv_3D);
  v.Register_Boundary(fd_bv_3D);
  w.Register_Boundary(fd_bv_3D);
  txx.Register_Boundary(fd_bv_3D);
  tyy.Register_Boundary(fd_bv_3D);
  tzz.Register_Boundary(fd_bv_3D);
  tyz.Register_Boundary(fd_bv_3D);
  txz.Register_Boundary(fd_bv_3D);
  txy.Register_Boundary(fd_bv_3D);

  fd_3D.Register_Array(u);
  fd_3D.Register_Array(v);
  fd_3D.Register_Array(w);
  fd_3D.Register_Array(txx);
  fd_3D.Register_Array(tyy);
  fd_3D.Register_Array(tzz);
  fd_3D.Register_Array(tyz);
  fd_3D.Register_Array(txz);
  fd_3D.Register_Array(txy);
  fd_3D.Register_Array(Lambda);
  fd_3D.Register_Array(Mu);
  fd_3D.Register_Array(Buoyancy);

  fd_3D.Register_Domain(I, J, K);
  u.Register_Shape(fd_shape_3D);
  v.Register_Shape(fd_shape_3D);
  w.Register_Shape(fd_shape_3D);

  
  // Initialization of prognostic fields
  for(int k=0;k<dimz;++k){
    for(int j=0;j<dimy;++j){
      for(int i=0;i<dimx;++i){
        u(0,k,j,i) = 0.0;
        v(0,k,j,i) = 0.0;
	w(0,k,j,i) = 0.0;

        txx(0,k,j,i) = 0.0;
        tyy(0,k,j,i) = 0.0;
        tzz(0,k,j,i) = 0.0;
	  
        tyz(0,k,j,i) = 0.0;
        txz(0,k,j,i) = 0.0;
        txy(0,k,j,i) = 0.0;
	  
        Buoyancy(0,k,j,i) = 1.0/((float (*)[dimy][dimx])rho.data())[k][j][i];
        Lambda(0,k,j,i) = ((float (*)[dimy][dimx])lam.data())[k][j][i];
        Mu(0,k,j,i) = ((float (*)[dimy][dimx])mu.data())[k][j][i];
      }
    }
  }
  
  
  Pochoir_Kernel_3D(fd_3D_velocity, t, k, j, i)
    // Update velocity
    u(t+1,k,j,i) = u(t,k,j,i) +
                   dt*Buoyancy(t,k,j,i)*(c0*(txx(t,k,j,i+1)-txx(t,k,j,i) + txy(t,k,j,i)-txy(t,k,j-1,i) + txz(t,k,j,i)-txz(t,k-1,j,i))
                   -c1*(txx(t,k,j,i+2)-txx(t,k,j,i-1) + txy(t,k,j+1,i)-txy(t,k,j-2,i) + txz(t,k+1,j,i)-txz(t,k-2,j,i)))/h;

    v(t+1,k,j,i) = v(t,k,j,i) +
                   dt*Buoyancy(t,k,j,i)*(c0*(txy(t,k,j,i)-txy(t,k,j,i-1) + tyy(t,k,j+1,i)-tyy(t,k,j,i) + tyz(t,k,j,i)-tyz(t,k-1,j,i))
                   -c1*(txy(t,k,j,i+1)-txy(t,k,j,i-2) + tyy(t,k,j+2,i)-tyy(t,k,j-1,i) + tyz(t,k+1,j,i)-tyz(t,k-2,j,i)))/h;

    w(t+1,k,j,i) = w(t,k,j,i) +
                   dt*Buoyancy(t,k,j,i)*(c0*(txz(t,k,j,i)-txz(t,k,j,i-1) + tyz(t,k,j,i)-tyz(t,k,j-1,i) + tzz(t,k+1,j,i)-tzz(t,k,j,i))
                   -c1*(txz(t,k,j,i+1)-txz(t,k,j,i-2) + tyz(t,k,j+1,i)-tyz(t,k,j-2,i) + tzz(t,k+2,j,i)-tzz(t,k-1,j,i)))/h;
           
  Pochoir_Kernel_End

  Pochoir_Kernel_3D(fd_3D_stress, t, k, j, i)
    // Update stress
    txx(t+1,k,j,i) = txx(t,k,j,i) + dt*((Lambda(t,k,j,i)+2.*Mu(t,k,j,i))*(c0*(u(t,k,j,i)-u(t,k,j,i-1))-c1*(u(t,k,j,i+1)-u(t,k,j,i-2)))
                                      + Lambda(t,k,j,i)*(c0*(v(t,k,j,i)-v(t,k,j-1,i) + w(t,k,j,i)-w(t,k-1,j,i)) - c1*(v(t,k,j+1,i)-v(t,k,j-2,i) + w(t,k+1,j,i)-w(t,k-2,j,i))))/h;
  
    tyy(t+1,k,j,i) = tyy(t,k,j,i) + dt*((Lambda(t,k,j,i)+2.*Mu(t,k,j,i))*(c0*(v(t,k,j,i)-v(t,k,j-1,i))-c1*(v(t,k,j+1,i)-v(t,k,j-2,i)))
                                      + Lambda(t,k,j,i)*(c0*(u(t,k,j,i)-u(t,k,j,i-1) + w(t,k,j,i)-w(t,k-1,j,i)) - c1*(u(t,k,j,i+1)-u(t,k,j,i-2) + w(t,k+1,j,i)-w(t,k-2,j,i))))/h;
  
    tzz(t+1,k,j,i) = tzz(t,k,j,i) + dt*((Lambda(t,k,j,i)+2.*Mu(t,k,j,i))*(c0*(w(t,k,j,i)-w(t,k-1,j,i))-c1*(w(t,k+1,j,i)-w(t,k-2,j,i)))
                                      + Lambda(t,k,j,i)*(c0*(u(t,k,j,i)-u(t,k,j,i-1) + v(t,k,j,i)-v(t,k,j-1,i)) - c1*(u(t,k,j,i+1)-u(t,k,j,i-2) + v(t,k,j+1,i)-v(t,k,j-2,i))))/h;
  
    tyz(t+1,k,j,i) = tyz(t,k,j,i) + dt*(Mu(t,k,j,i)*(c0*(v(t,k+1,j,i)-v(t,k,j,i) + w(t,k,j+1,i)-w(t,k,j,i)) - c1*(v(t,k+2,j,i)-v(t,k-1,j,i) + w(t,k,j+2,i)-w(t,k,j-1,i))))/h;
  
    txz(t+1,k,j,i) = txz(t,k,j,i) + dt*(Mu(t,k,j,i)*(c0*(u(t,k+1,j,i)-u(t,k,j,i) + w(t,k,j,i+1)-w(t,k,j,i)) - c1*(u(t,k+2,j,i)-u(t,k-1,j,i) + w(t,k,j,i+2)-w(t,k,j,i-1))))/h;
  
    txy(t+1,k,j,i) = txy(t,k,j,i) + dt*(Mu(t,k,j,i)*(c0*(u(t,k,j+1,i)-u(t,k,j,i) + v(t,k,j,i+1)-v(t,k,j,i)) - c1*(u(t,k,j+2,i)-u(t,k,j-1,i) + v(t,k,j,i+2)-v(t,k,j,i-1))))/h;
    
  Pochoir_Kernel_End
  
  Pochoir_Kernel_3D(fd_3D_velocity_swap, t, k, j, i)
    // swap velocity result   
    u(t,k,j,i) = u(t+1,k,j,i);                   
    v(t,k,j,i) = v(t+1,k,j,i);
    w(t,k,j,i) = w(t+1,k,j,i);
    
  Pochoir_Kernel_End

  Pochoir_Kernel_3D(fd_3D_stress_swap, t, k, j, i)
    // swap stress result
    txx(t,k,j,i) = txx(t+1,k,j,i);
    tyy(t,k,j,i) = tyy(t+1,k,j,i);
    tzz(t,k,j,i) = tzz(t+1,k,j,i);
    tyz(t,k,j,i) = tyz(t+1,k,j,i);
    txz(t,k,j,i) = txz(t+1,k,j,i);
    txy(t,k,j,i) = txy(t+1,k,j,i);
    
  Pochoir_Kernel_End

  // Location of source.
  int sx=(int)round(coorsrc[0]/h);
  int sy=(int)round(coorsrc[1]/h);
  int sz=(int)round(coorsrc[2]/h);
   // Set up solution fields.
  std::vector<float> u_ref(dimx*dimy*dimz), v_ref(dimx*dimy*dimz), w_ref(dimx*dimy*dimz),
    txx_ref(dimx*dimy*dimz), tyy_ref(dimx*dimy*dimz), tzz_ref(dimx*dimy*dimz);
    //,tyz_ref(dimx*dimy*dimz), txz_ref(dimx*dimy*dimz), txy_ref(dimx*dimy*dimz);
    
  std::vector<float> uss, vss, wss, pss;
  uss.reserve(nrec*ntsteps);
  vss.reserve(nrec*ntsteps);
  wss.reserve(nrec*ntsteps);
  pss.reserve(nrec*ntsteps);

  

  for(int times=0; times<2*ntsteps;++times){
    if(times%2==0){
    
      	fd_3D.Run(1, fd_3D_stress);      
        fd_3D.Run(1, fd_3D_stress_swap);      
     
    }else{
     
        fd_3D.Run(1, fd_3D_velocity);      
      	fd_3D.Run(1, fd_3D_velocity_swap);
      	
  #pragma omp parallel
  {       	
      #pragma omp single nowait
      {
	for(int i=0;i<nrec;i++){
	  int xi = (int)round(coorrec[i*3]/h);
	  int yi = (int)round(coorrec[i*3+1]/h);
	  int zi = (int)round(coorrec[i*3+2]/h);

	 uss.push_back(u(1,zi,yi,xi));
	}
      }

      #pragma omp single nowait
      {
	for(int i=0;i<nrec;i++){
	  int xi = (int)round(coorrec[i*3]/h);
	  int yi = (int)round(coorrec[i*3+1]/h);
	  int zi = (int)round(coorrec[i*3+2]/h);
	  
	  vss.push_back(v(1,zi,yi,xi));
	}
      }

      #pragma omp single nowait
      {
	for(int i=0;i<nrec;i++){
	  int xi = (int)round(coorrec[i*3]/h);
	  int yi = (int)round(coorrec[i*3+1]/h);
	  int zi = (int)round(coorrec[i*3+2]/h);
	  
	  wss.push_back(w(1,zi,yi,xi));
	}
      }
      
      #pragma omp single
      {
	for(int i=0;i<nrec;i++){
	  int xi = (int)round(coorrec[i*3]/h);
	  int yi = (int)round(coorrec[i*3+1]/h);
	  int zi = (int)round(coorrec[i*3+2]/h);

	  pss.push_back((txx(1,zi,yi,xi)+
			 tyy(1,zi,yi,xi)+
			 tzz(1,zi,yi,xi))/3);
	}
      }
      
      #pragma omp single
      {
        if(times/2<snt){ // Add source        
	   txx(0,sz,sy,sx) = txx(1,sz,sy,sx) -= xsrc[times/2]/3;
	   tyy(0,sz,sy,sx) = tyy(1,sz,sy,sx) -= ysrc[times/2]/3;
	   tzz(0,sz,sy,sx) = tzz(1,sz,sy,sx) -= zsrc[times/2]/3;
        }
      }    
      
     }//end of parallel region          
    }
  }
  
  for(int k=0;k<dimz;++k){
    for(int j=0;j<dimy;++j){
      for(int i=0;i<dimx;++i){
       
        index = k*dimx*dimy+j*dimz+i;
        
        u_ref[index] = u.interior(1,k,j,i);
        v_ref[index] = v.interior(1,k,j,i);
        w_ref[index] = w.interior(1,k,j,i);
        
        txx_ref[index] = txx.interior(1,k,j,i);
        tyy_ref[index] = tyy.interior(1,k,j,i);
        tzz_ref[index] = tzz.interior(1,k,j,i);
        //tyz_ref[index] = tyz.interior(1,k,j,i);
        //txz_ref[index] = txz.interior(1,k,j,i);
        //txy_ref[index] = txy.interior(1,k,j,i);
         
    }
   }
  }
  
  {
    int dims[]={dimx, dimy, dimz};
    float spacing[]={h, h, h};
    opesci_dump_solution_vts("solution_pochoir", dims, spacing, u_ref,v_ref,w_ref,txx_ref,tyy_ref,tzz_ref);
  }
  
  {
    int dims[]={(int)round(sqrt(nrec)), (int)round(sqrt(nrec)), ntsteps};
    float spacing[]={h, h, dt};
    opesci_dump_receivers_vts("receivers_pochoir", dims, spacing, uss, vss, wss, pss);
  }

  return 0;
}
