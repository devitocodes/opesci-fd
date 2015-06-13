<%include file="copyright.txt"/>

<%include file="common_include.txt"/>

int main(){

<%include file="read_data.txt"/>

  // time periodicity for update
  const int _tp = ${time_period};

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

    float (*U)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _u.data();
    float (*V)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _v.data();
    float (*W)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _w.data();

    float (*Txx)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _txx.data();
    float (*Tyy)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _tyy.data();
    float (*Tzz)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _tzz.data();
    float (*Tyz)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _tyz.data();
    float (*Txz)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _txz.data();
    float (*Txy)[dimy][dimz][_tp] = (float (*)[dimy][dimz][_tp]) _txy.data();

    // main time loop
    for(int _ti=0;_ti<ntsteps;_ti++){

      // shared variables
      int t = _ti % _tp; // array index of current time step
      int t1 = (t+1) % _tp; // array index of the grid to be updated

      // Compute stresses
    #pragma omp for schedule(guided)
      for(int x=2;x<dimx-2;x++){
      	for(int y=2;y<dimy-2;y++){
      	  for(int z=2;z<dimz-2;z++){
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
        for(int y=2;y<dimy-2;y++){
          for(int z=2;z<dimz-2;z++){
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

      	  pss.push_back((Txx[xi][yi][zi][t1] + Tyy[xi][yi][zi][t1] + Tzz[xi][yi][zi][t1])/3);
      	}
      }

      // Insert source
    #pragma omp single
      {
      	if(_ti<snt){
      	  Txx[_sx][_sy][_sz][t1] -= xsrc[_ti]/3;
      	  Tyy[_sx][_sy][_sy][t1] -= ysrc[_ti]/3;
      	  Tzz[_sx][_sy][_sy][t1] -= zsrc[_ti]/3;
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

  {
    int dims[]={dimx, dimy, dimz};
    float spacing[]={h, h, h};
    opesci_dump_solution_vts("solution_tmpl", dims, spacing, _u_out,_v_out,_w_out,_txx_out,_tyy_out,_tzz_out);
  }

  // output receiver reading
  {
    int dims[]={(int)round(sqrt(nrec)), (int)round(sqrt(nrec)), ntsteps};
    float spacing[]={h, h, dt};
    opesci_dump_receivers_vts("receivers_tmpl", dims, spacing, uss, vss, wss, pss);    
  }

  return 0;
}
