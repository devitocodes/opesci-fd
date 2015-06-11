// insert copyright information

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
  std::vector<float> _rho;
  if(opesci_read_simple_binary(rhofile.c_str(), _rho)){
    opesci_abort("Missing input file\n");
  }
  assert(rho.size()==dimx*dimy*dimz);
  std::vector<float> _buoyancy(rho.size());

  // Calculate Lame constents.
  std::vector<float> _mu, _lam;
  opesci_calculate_lame_costants(vp, vs, rho, _mu, _lam);

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

  // time periodicity for update
  int _tp = ${time_period}

  // Set up solution fields.
  std::vector<float> _u(dimx*dimy*dimz*_tp), _v(dimx*dimy*dimz*_tp), _w(dimx*dimy*dimz*_tp),
    _txx(dimx*dimy*dimz*_tp), _tyy(dimx*dimy*dimz*_tp), _tzz(dimx*dimy*dimz*_tp),
    _tyz(dimx*dimy*dimz*_tp), _txz(dimx*dimy*dimz*_tp), _txy(dimx*dimy*dimz*_tp);

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
      _u[i*_tp] = 0.0;
      _v[i*_tp] = 0.0;
      _w[i*_tp] = 0.0;

      _txx[i*_tp] = 0.0;
      _tyy[i*_tp] = 0.0;
      _tzz[i*_tp] = 0.0;
      _tyz[i*_tp] = 0.0;
      _txz[i*_tp] = 0.0;
      _txy[i*_tp] = 0.0;

      _buoyancy[i] = 1.0/_rho[i];
    }

    // initialise at t=0
    float (*lambda)[dimy][dimz] = (float (*)[dimy][dimz]) _lam.data();
    float (*mu)[dimy][dimz] = (float (*)[dimy][dimz]) _mu.data(); // need to correct to use the effective media parameters
    float (*beta)[dimy][dimz] = (float (*)[dimy][dimz]) _buoyancy.data(); // need to correct to use the effective media parameters

    float (*U)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _u.data();
    float (*V)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _v.data();
    float (*W)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _w.data();

    float (*Txx)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _txx.data();
    float (*Tyy)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _tyy.data();
    float (*Tzz)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _tzz.data();
    float (*Tyz)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _tyz.data();
    float (*Txz)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _txz.data();
    float (*Txy)[dimy][dimz][_tp] = (float (*)[dimy][dimz]) _txy.data();

    // main time loop
    for(int _ti=0;_ti<ntsteps;_ti++){

      // shared variables
    #pragma omp single
      {
        int t = _ti % _tp; // array index of current time step
        int t1 = (t+1) % _tp // array index of the grid to be updated
      }

      // Compute stresses
    #pragma omp for schedule(guided)
      for(int x=2;x<dimx-2;x++){
      	for(int y=2;j<dimy-2;y++){
      	  for(int z=2;i<dimz-2;z++){
      	    ${Txx};
      	    ${Tyy};
      	    ${Tzz};
      	    ${Tyz};
      	    ${Txz};
      	    ${Txy};
      	  }
      	}
      }

      // Compute velocities
    #pragma omp for schedule(guided)
      for(int x=2;x<dimx-2;x++){
        for(int y=2;j<dimy-2;y++){
          for(int z=2;i<dimz-2;z++){
            ${U};
            ${V};
            ${W};
          }
        }
      }
   
#pragma omp single nowait
      {
      	for(int i=0;i<nrec;i++){
      	  int xi = (int)round(coorrec[i*3]/h);
      	  int yi = (int)round(coorrec[i*3+1]/h);
      	  int zi = (int)round(coorrec[i*3+2]/h);

      	  uss.push_back(U[xi][yi][zi][t1]);
      	}
      }

#pragma omp single nowait
      {
      	for(int i=0;i<nrec;i++){
      	  int xi = (int)round(coorrec[i*3]/h);
      	  int yi = (int)round(coorrec[i*3+1]/h);
      	  int zi = (int)round(coorrec[i*3+2]/h);
      	  
      	  vss.push_back(V[xi][yi][zi][t1]);
      	}
      }

#pragma omp single nowait
      {
      	for(int i=0;i<nrec;i++){
      	  int xi = (int)round(coorrec[i*3]/h);
      	  int yi = (int)round(coorrec[i*3+1]/h);
      	  int zi = (int)round(coorrec[i*3+2]/h);
      	  
      	  wss.push_back(W[xi][yi][zi][t1]);
      	}
      }

      // A barrier is implied here because the next block modifies the stress field.
#pragma omp single
      {
      	for(int i=0;i<nrec;i++){
      	  int xi = (int)round(coorrec[i*3]/h);
      	  int yi = (int)round(coorrec[i*3+1]/h);
      	  int zi = (int)round(coorrec[i*3+2]/h);

      	  pss.push_back((Txx[zi][yi][xi][t1] + Tyy[zi][yi][xi][t1] + Tzz[zi][yi][xi][t1])/3);
      	}
      }

      // Insert source
#pragma omp single
      {
      	if(ti<snt){
      	  int _sx=(int)round(coorsrc[0]/h);
      	  int _sy=(int)round(coorsrc[1]/h);
      	  int _sz=(int)round(coorsrc[2]/h);

      	  assert(_sx<dimx);
      	  assert(_sy<dimy);
      	  assert(_sz<dimz);

      	  Txx[_sx][_sy][_sz][t1] -= xsrc[_ti]/3;
      	  Tyy[_sz][_sy][_sx][t1] -= ysrc[_ti]/3;
      	  Tzz[_sz][_sy][_sx][t1] -= zsrc[_ti]/3;
      	}
      }
    } // end of time loop
  } // end of parallel section

  // output solution
  std::vector<float> _u_out(dimx*dimy*dimz), _v_out(dimx*dimy*dimz), _w_out(dimx*dimy*dimz),
    _txx_out(dimx*dimy*dimz), _tyy_out(dimx*dimy*dimz), _tzz_out(dimx*dimy*dimz);

  int _last_loop = ntsteps % _tp;
  for(int i=0;i<dimx*dimy*dimz;i++){
    _u_out[i] = _u[i*_tp+_last_loop];
    _v_out[i] = _v[i*_tp+_last_loop];
    _w_out[i] = _w[i*_tp+_last_loop];
    _txx_out[i] = _txx[i*_tp+_last_loop];
    _tyy_out[i] = _tyy[i*_tp+_last_loop];
    _tzz_out[i] = _tzz[i*_tp+_last_loop];
  }

  int dims[]={dimx, dimy, dimz};
  float spacing[]={h, h, h};
  opesci_dump_solution_vts("solution_tmpl", dims, spacing, _u_out,_v_out,_w_out,_txx_out,_tyy_out,_tzz_out);
  
  // output receiver reading
  int dims[]={(int)round(sqrt(nrec)), (int)round(sqrt(nrec)), ntsteps};
  float spacing[]={h, h, dt};
  opesci_dump_receivers_vts("receivers_tmpl", dims, spacing, uss, vss, wss, pss);

  return 0;
}
