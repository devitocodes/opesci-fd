<%include file="copyright.txt"/>

<%include file="common_include.txt"/>
#include <cmath>

int main(){

  float rho = 1.0;
  float mu = 0.5;
  float lambda = 0.5;
  float Vp = sqrt((lamda+2*mu)/rho);
  float Vs = sqrt(mu/rho);
  float h = 0.025;
  float dt = 0.5*h/Vp/2; // need to check
  int margin = 2; // ghost cells for boundary conditions
  int dimx, dimy;
  dimx = dimy = (int) 1.0/h + margin*2;

  // time periodicity for update
  const int _tp = ${time_period};

  // set up solution mesh
  float U[_tp][dimx][dimy];
  float V[_tp][dimx][dimy];
  float Txx[_tp][dimx][dimy];
  float Tyy[_tp][dimx][dimy];
  float Txy[_tp][dimx][dimy];

  // Set up seismic sections
  std::vector<float> uss, vss, wss, pss;
  uss.reserve(nrec*ntsteps);
  vss.reserve(nrec*ntsteps);
  wss.reserve(nrec*ntsteps);
  pss.reserve(nrec*ntsteps);

  // shared variables
  int t, t1;

#pragma omp parallel
  {
    // Initialise fields exploiting first touch.
    // main time loop
  #pragma omp for
  for(int _ti=0;_ti<ntsteps;_ti++){

      // shared variables
      int t = _ti % _tp; // array index of current time step
      int t1 = (t+1) % _tp; // array index of the grid to be updated

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
      Tyy[t1][x][dimy+3] = -Tyy[t1][x][dimy+1];
      Txy[t1][x][dimy+3] = -Txy[t1][x][dimy];
      Txy[t1][x][dimy+2] = -Txy[t1][x][dimy+1];
    }
    #pragma omp for
    for(int y=0;y<dimy;y++){
      // boundary x=2
      Tyy[t1][1][y] = -Tyy[t1][3][y];
      Txy[t1][0][y] = -Txy[t1][3][y];
      Txy[t1][1][y] = -Txy[t1][2][y];
      // boundary x=dimx+2
      Tyy[t1][dimx+3][y] = -Tyy[t1][dimx+1][y];
      Txy[t1][dimx+3][y] = -Txy[t1][dimx][y];
      Txy[t1][dimx+2][y] = -Txy[t1][dimx+1][y];
    }

    // Compute velocities
    #pragma omp for
    for(int x=2;x<dimx-2;x++){
      for(int y=2;y<dimy-2;y++){
        ${U};
        ${V};
      }
    }

    // update ghost cells for boundary conditions
    #pragma omp for
    for(int x=margin;x<dimx-margin;x++){
      for(int y=margin;y<dimy-margin;y++){
        ${Txx};
        ${Tyy};
        ${Txy};
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
