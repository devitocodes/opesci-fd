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

int main(){
  std::string vpfile("VPhomogx200");       // Vp file in binary (in m/s)
  std::string vsfile("VShomogx200");       // Vs file: Vs file in binary (in m/s)
  std::string rhofile("RHOhomogx200");     // rhofile: densities file in binary (in kg/m**3)
  std::string geomsrc("geomx200.src");     // Geometry file with sources locations (in m)
  std::string geomrec("geomx200.rec");     // Geometry file with receivers locations (in m)
  int dimx=200, dimy=200, dimz=200;   // Model dimensions in number of nodes
  float h=25;                         // Node spacing in meters (hx=hy=hz=h)
  std::string xsrcname("xrick5.sep");      // Source files
  std::string ysrcname("yrick5.sep");
  std::string zsrcname("zrick5.sep"); 
  float sdt=0.004;                    // Source sample rate (in s).
  float maxt=1.0;                     // Maximum time to compute;
  std::string useis("seisu");              // Seismogram files for the x, y and z components of velocity.
  std::string vseis("seisv");
  std::string wseis("seisw");
  std::string pseis("seisp");              // Seismogram files for pressure.
  
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

  // Set up solution fields.
  std::vector<float> u(dimx*dimy*dimz), v(dimx*dimy*dimz), w(dimx*dimy*dimz),
    txx(dimx*dimy*dimz), tyy(dimx*dimy*dimz), tzz(dimx*dimy*dimz),
    tyz(dimx*dimy*dimz), txz(dimx*dimy*dimz), txy(dimx*dimy*dimz);

  // Set up seismic sections
  std::vector<float> uss, vss, wss, pss;
  uss.reserve(nrec*ntsteps);
  vss.reserve(nrec*ntsteps);
  wss.reserve(nrec*ntsteps);
  pss.reserve(nrec*ntsteps);

#pragma omp parallel
  {
    // Initialise fields exploiting first touch.
#pragma omp for nowait
    for(int i=0;i<dimx*dimy*dimz;i++){
      u[i] = 0.0;
      v[i] = 0.0;
      w[i] = 0.0;

      txx[i] = 0.0;
      tyy[i] = 0.0;
      tzz[i] = 0.0;
      tyz[i] = 0.0;
      txz[i] = 0.0;
      txy[i] = 0.0;

      buoyancy[i] = 1.0/rho[i];
    }

    // Start of time loop.
    float (*Lambda)[dimy][dimx] = (float (*)[dimy][dimx]) lam.data();
    float (*Mu)[dimy][dimx] = (float (*)[dimy][dimx]) mu.data(); // need to correct to use the effective media parameters
    float (*Buoyancy)[dimy][dimx] = (float (*)[dimy][dimx]) buoyancy.data(); // need to correct to use the effective media parameters

    float (*U)[dimy][dimx] = (float (*)[dimy][dimx]) u.data();
    float (*V)[dimy][dimx] = (float (*)[dimy][dimx]) v.data();
    float (*W)[dimy][dimx] = (float (*)[dimy][dimx]) w.data();

    float (*Txx)[dimy][dimx] = (float (*)[dimy][dimx]) txx.data();
    float (*Tyy)[dimy][dimx] = (float (*)[dimy][dimx]) tyy.data();
    float (*Tzz)[dimy][dimx] = (float (*)[dimy][dimx]) tzz.data();
    float (*Tyz)[dimy][dimx] = (float (*)[dimy][dimx]) tyz.data();
    float (*Txz)[dimy][dimx] = (float (*)[dimy][dimx]) txz.data();
    float (*Txy)[dimy][dimx] = (float (*)[dimy][dimx]) txy.data();

    float c0=9./8.;
    float c1=1./24.;

    for(int ti=0;ti<ntsteps;ti++){
      
      // Compute stresses
#pragma omp for
      for(int k=2;k<dimz-2;k++){
	for(int j=2;j<dimy-2;j++){
	  for(int i=2;i<dimx-2;i++){
	    Txx[k][j][i] += dt*((Lambda[k][j][i]+2.*Mu[k][j][i])*(c0*(U[k][j][i]-U[k][j][i-1])-c1*(U[k][j][i+1]-U[k][j][i-2]))
				+ Lambda[k][j][i]*(c0*(V[k][j][i]-V[k][j-1][i] + W[k][j][i]-W[k-1][j][i]) - c1*(V[k][j+1][i]-V[k][j-2][i] + W[k+1][j][i]-W[k-2][j][i])))/h;
	    
	    Tyy[k][j][i] += dt*((Lambda[k][j][i]+2.*Mu[k][j][i])*(c0*(V[k][j][i]-V[k][j-1][i])-c1*(V[k][j+1][i]-V[k][j-2][i]))
				+ Lambda[k][j][i]*(c0*(U[k][j][i]-U[k][j][i-1] + W[k][j][i]-W[k-1][j][i]) - c1*(U[k][j][i+1]-U[k][j][i-2] + W[k+1][j][i]-W[k-2][j][i])))/h;
	    
	    Tzz[k][j][i] += dt*((Lambda[k][j][i]+2.*Mu[k][j][i])*(c0*(W[k][j][i]-W[k-1][j][i])-c1*(W[k+1][j][i]-W[k-2][j][i]))
				+ Lambda[k][j][i]*(c0*(U[k][j][i]-U[k][j][i-1] + V[k][j][i]-V[k][j-1][i]) - c1*(U[k][j][i+1]-U[k][j][i-2] + V[k][j+1][i]-V[k][j-2][i])))/h;
	    
	    Tyz[k][j][i] += dt*(Mu[k][j][i]*(c0*(V[k+1][j][i]-V[k][j][i] + W[k][j+1][i]-W[k][j][i]) - c1*(V[k+2][j][i]-V[k-1][j][i] + W[k][j+2][i]-W[k][j-1][i])))/h;
	    
	    Txz[k][j][i] += dt*(Mu[k][j][i]*(c0*(U[k+1][j][i]-U[k][j][i] + W[k][j][i+1]-W[k][j][i]) - c1*(U[k+2][j][i]-U[k-1][j][i] + W[k][j][i+2]-W[k][j][i-1])))/h;
	    
	    Txy[k][j][i] += dt*(Mu[k][j][i]*(c0*(U[k][j+1][i]-U[k][j][i] + V[k][j][i+1]-V[k][j][i]) - c1*(U[k][j+2][i]-U[k][j-1][i] + V[k][j][i+2]-V[k][j][i-1])))/h;

	    if(!std::isfinite(Tyz[k][j][i])){
	      
	      std::cerr<<"a="<<(c0*(V[k+1][j][i]-V[k][j][i] + W[k][j+1][i]-W[k][j][i])
				-c1*(V[k+2][j][i]-V[k-1][j][i] + W[k][j+2][i]-W[k][j-1][i]))<<std::endl
		       <<V[k+1][j][i]<<", "<<V[k][j][i]<<", "<<W[k][j+1][i]<<", "<<W[k][j][i]<<std::endl
		       <<V[k+2][j][i]<<", "<<V[k-1][j][i]<<", "<<W[k][j+2][i]<<", "<<W[k][j-1][i]<<std::endl;
	      std::cerr<<"b="<<dt*Mu[k][j][i]<<std::endl;
	    }
	  }
	}
      }

      // Compute velocities
#pragma omp for
      for(int k=2;k<dimz-2;k++){
	for(int j=2;j<dimy-2;j++){
	  for(int i=2;i<dimx-2;i++){
	    
	    U[k][j][i] +=
	      dt*Buoyancy[k][j][i]*(c0*(Txx[k][j][i+1]-Txx[k][j][i] + Txy[k][j][i]-Txy[k][j-1][i] + Txz[k][j][i]-Txz[k-1][j][i])
				    -c1*(Txx[k][j][i+2]-Txx[k][j][i-1] + Txy[k][j+1][i]-Txy[k][j-2][i] + Txz[k+1][j][i]-Txz[k-2][j][i]))/h;
	    
	    V[k][j][i] +=
	      dt*Buoyancy[k][j][i]*(c0*(Txy[k][j][i]-Txy[k][j][i-1] + Tyy[k][j+1][i]-Tyy[k][j][i] + Tyz[k][j][i]-Tyz[k-1][j][i])
				    -c1*(Txy[k][j][i+1]-Txy[k][j][i-2] + Tyy[k][j+2][i]-Tyy[k][j-1][i] + Tyz[k+1][j][i]-Tyz[k-2][j][i]))/h;
	    
	    W[k][j][i] +=
	      dt*Buoyancy[k][j][i]*(c0*(Txz[k][j][i]-Txz[k][j][i-1] + Tyz[k][j][i]-Tyz[k][j-1][i] + Tzz[k+1][j][i]-Tzz[k][j][i])
				    -c1*(Txz[k][j][i+1]-Txz[k][j][i-2] + Tyz[k][j+1][i]-Tyz[k][j-2][i] + Tzz[k+2][j][i]-Tzz[k-1][j][i]))/h;
	  }
	}
      }
   
#pragma omp single nowait
      {
	for(int i=0;i<nrec;i++){
	  int xi = (int)round(coorrec[i*3]/h);
	  int yi = (int)round(coorrec[i*3+1]/h);
	  int zi = (int)round(coorrec[i*3+2]/h);

	  uss.push_back(U[zi][yi][xi]);
	}
      }

#pragma omp single nowait
      {
	for(int i=0;i<nrec;i++){
	  int xi = (int)round(coorrec[i*3]/h);
	  int yi = (int)round(coorrec[i*3+1]/h);
	  int zi = (int)round(coorrec[i*3+2]/h);
	  
	  vss.push_back(V[zi][yi][xi]);
	}
      }

#pragma omp single nowait
      {
	for(int i=0;i<nrec;i++){
	  int xi = (int)round(coorrec[i*3]/h);
	  int yi = (int)round(coorrec[i*3+1]/h);
	  int zi = (int)round(coorrec[i*3+2]/h);
	  
	  wss.push_back(W[zi][yi][xi]);
	}
      }

      // A barrier is implied here because the next block modifies the stress field.
#pragma omp single
      {
	for(int i=0;i<nrec;i++){
	  int xi = (int)round(coorrec[i*3]/h);
	  int yi = (int)round(coorrec[i*3+1]/h);
	  int zi = (int)round(coorrec[i*3+2]/h);

	  pss.push_back((Txx[zi][yi][xi]+
			 Tyy[zi][yi][xi]+
			 Tzz[zi][yi][xi])/3);
	}
      }

      // Insert source
#pragma omp single
      {
	if(ti<snt){
	  int sx=(int)round(coorsrc[0]/h);
	  int sy=(int)round(coorsrc[1]/h);
	  int sz=(int)round(coorsrc[2]/h);

	  assert(sx<dimx);
	  assert(sy<dimy);
	  assert(sz<dimz);

	  Txx[sz][sy][sx] -= xsrc[ti]/3;
	  Tyy[sz][sy][sx] -= ysrc[ti]/3;
	  Tzz[sz][sy][sx] -= zsrc[ti]/3;
	}
      }
    }
  }

  {
    int dims[]={dimx, dimy, dimz};
    float spacing[]={h, h, h};
    opesci_dump_solution_vts("solution", dims, spacing, u,v,w,txx,tyy,tzz);
  }
  
  {
    int dims[]={(int)round(sqrt(nrec)), (int)round(sqrt(nrec)), ntsteps};
    float spacing[]={h, h, dt};
    opesci_dump_field_vts("receivers", dims, spacing, uss);
  }

  std::ofstream ussfh;
  ussfh.open ("useis.raw", std::ios::out | std::ios::trunc | std::ios::binary); 
  ussfh.write((char *)uss.data(), uss.size()*sizeof(float));
  ussfh.close();

  std::ofstream vssfh;
  vssfh.open ("vseis.raw", std::ios::out | std::ios::trunc | std::ios::binary); 
  vssfh.write((char *)vss.data(), vss.size()*sizeof(float));
  vssfh.close();

  std::ofstream wssfh;
  wssfh.open ("wseis.raw", std::ios::out | std::ios::trunc | std::ios::binary); 
  wssfh.write((char *)wss.data(), wss.size()*sizeof(float));
  wssfh.close();

  std::ofstream pssfh;
  pssfh.open ("pseis.raw", std::ios::out | std::ios::trunc | std::ios::binary); 
  pssfh.write((char *)pss.data(), pss.size()*sizeof(float));
  pssfh.close();

  return 0;
}
