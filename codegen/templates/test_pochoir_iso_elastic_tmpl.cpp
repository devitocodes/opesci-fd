<%include file="copyright.txt"/>

<%include file="common_include.txt"/>

#include <pochoir.hpp>

Pochoir_Boundary_3D(fd_bv_3D, arr, t, x, y, z)
    return 0;
Pochoir_Boundary_End

int main(){

<%include file="read_data.txt"/>

  // time periodicity for update
  const int _tp = ${time_period};
  
  Pochoir_Shape_3D fd_shape_3D[] = {
    {1,0,0,0},
    {0,0,0,0}
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
  Pochoir_Array<float, 3> U(dimx, dimy, dimz), V(dimx, dimy, dimz), W(dimx, dimy, dimz),
    Txx(dimx, dimy, dimz), Tyy(dimx, dimy, dimz), Tzz(dimx, dimy, dimz),
    Tyz(dimx, dimy, dimz), Txz(dimx, dimy, dimz), Txy(dimx, dimy, dimz);

  // Subsurface model
  Pochoir_Array<float, 3> beta(dimx, dimy, dimz), lambda(dimx, dimy, dimz), mu(dimx, dimy, dimz);
  
  Pochoir_3D fd_3D(fd_shape_3D);
  int ds=2;
  Pochoir_Domain I(0+ds, dimx-ds), J(0+ds, dimy-ds), K(0+ds, dimz-ds);
  
  U.Register_Boundary(fd_bv_3D);
  V.Register_Boundary(fd_bv_3D);
  W.Register_Boundary(fd_bv_3D);
  Txx.Register_Boundary(fd_bv_3D);
  Tyy.Register_Boundary(fd_bv_3D);
  Tzz.Register_Boundary(fd_bv_3D);
  Tyz.Register_Boundary(fd_bv_3D);
  Txz.Register_Boundary(fd_bv_3D);
  Txy.Register_Boundary(fd_bv_3D);

  fd_3D.Register_Array(U);
  fd_3D.Register_Array(V);
  fd_3D.Register_Array(W);
  fd_3D.Register_Array(Txx);
  fd_3D.Register_Array(Tyy);
  fd_3D.Register_Array(Tzz);
  fd_3D.Register_Array(Tyz);
  fd_3D.Register_Array(Txz);
  fd_3D.Register_Array(Txy);
  fd_3D.Register_Array(lambda);
  fd_3D.Register_Array(mu);
  fd_3D.Register_Array(beta);

  fd_3D.Register_Domain(I, J, K);
  U.Register_Shape(fd_shape_3D);
  V.Register_Shape(fd_shape_3D);
  W.Register_Shape(fd_shape_3D);
  
  // Initialization of prognostic fields
  for(int x=0;x<dimx;++x){
    for(int y=0;y<dimy;++y){
      for(int z=0;z<dimx;++z){
        U(0,x,y,z) = 0.0;
        V(0,x,y,z) = 0.0;
        W(0,x,y,z) = 0.0;

        Txx(0,x,y,z) = 0.0;
        Tyy(0,x,y,z) = 0.0;
        Tzz(0,x,y,z) = 0.0;

        Tyz(0,x,y,z) = 0.0;
        Txz(0,x,y,z) = 0.0;
        Txy(0,x,y,z) = 0.0;
	  
        beta(0,x,y,z) = 1.0/((float (*)[dimy][dimz])rho.data())[x][y][z];
        lambda(0,x,y,z) = ((float (*)[dimy][dimz])lam.data())[x][y][z];
        mu(0,x,y,z) = ((float (*)[dimy][dimz])mu.data())[x][y][z];
      }
    }
  }
  
  
  Pochoir_Kernel_3D(fd_3D_velocity, t, x, y, z)
    // Update velocity
    ${U}
    ${V}
    ${W}
           
  Pochoir_Kernel_End

  Pochoir_Kernel_3D(fd_3D_stress, t, k, j, i)
    // Update stress
    ${Txx}
    ${Tyy}
    ${Tzz}
    ${Tyz}
    ${Txz}
    ${Txy}

  Pochoir_Kernel_End
  
  Pochoir_Kernel_3D(fd_3D_velocity_swap, t, x, y, z)
    // swap velocity result   
    U(t,x,y,z) = U(t+1,x,y,z);                   
    V(t,x,y,z) = V(t+1,x,y,z);
    W(t,x,y,z) = W(t+1,x,y,z);
    
  Pochoir_Kernel_End

  Pochoir_Kernel_3D(fd_3D_stress_swap, t, x, y, z)
    // swap stress result
    Txx(t,x,y,z) = Txx(t+1,x,y,z);
    Tyy(t,x,y,z) = Tyy(t+1,x,y,z);
    Tzz(t,x,y,z) = Tzz(t+1,x,y,z);
    Tyz(t,x,y,z) = Tyz(t+1,x,y,z);
    Txz(t,x,y,z) = Txz(t+1,x,y,z);
    Txy(t,x,y,z) = Txy(t+1,x,y,z);
    
  Pochoir_Kernel_End

  // Location of source.
  int sx=(int)round(coorsrc[0]/h);
  int sy=(int)round(coorsrc[1]/h);
  int sz=(int)round(coorsrc[2]/h);
   // Set up solution fields.
  std::vector<float> _u_out(dimx*dimy*dimz), _v_out(dimx*dimy*dimz), _w_out(dimx*dimy*dimz),
    _txx_out(dimx*dimy*dimz), _tyy_out(dimx*dimy*dimz), _tzz_out(dimx*dimy*dimz);
    
  std::vector<float> uss, vss, wss, pss;
  uss.reserve(nrec*ntsteps);
  vss.reserve(nrec*ntsteps);
  wss.reserve(nrec*ntsteps);
  pss.reserve(nrec*ntsteps);

  

  for(int _t=0; _t<2*ntsteps;++_t){
    if(_t%2==0){
    
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

          uss.push_back(U(1,xi,yi,zi));
        }
        }

        #pragma omp single nowait
        {
      	for(int i=0;i<nrec;i++){
      	  int xi = (int)round(coorrec[i*3]/h);
      	  int yi = (int)round(coorrec[i*3+1]/h);
      	  int zi = (int)round(coorrec[i*3+2]/h);

          vss.push_back(V(1,xi,yi,zi));
        }
        }

        #pragma omp single nowait
        {
      	for(int i=0;i<nrec;i++){
      	  int xi = (int)round(coorrec[i*3]/h);
      	  int yi = (int)round(coorrec[i*3+1]/h);
      	  int zi = (int)round(coorrec[i*3+2]/h);
      	  
      	  wss.push_back(W(1,xi,yi,zi));
        }
        }
        
        #pragma omp single
        {
      	for(int i=0;i<nrec;i++){
      	  int xi = (int)round(coorrec[i*3]/h);
      	  int yi = (int)round(coorrec[i*3+1]/h);
      	  int zi = (int)round(coorrec[i*3+2]/h);

      	  pss.push_back((Txx(1,xi,yi,zi)+Tyy(1,xi,yi,zi)+Tzz(1,xi,yi,zi))/3);
        }
        }
        
        // Add source
        #pragma omp single
        {
          if(_t/2<snt){
            Txx(0,sx,sy,sz) = Txx(1,sx,sy,sz) -= xsrc[_t/2]/3;
            Tyy(0,sx,sy,sz) = Tyy(1,sx,sy,sz) -= ysrc[_t/2]/3;
            Tzz(0,sx,sy,sz) = Tzz(1,sx,sy,sz) -= zsrc[_t/2]/3;
          }
        }
      } //end of parallel region
    } //end else
  } //end for
  
  for(int x=0;x<dimx;++x){
    for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
       
        index = x*dimy*dimz+y*dimz+z;
        
        _u_out[index] = U.interior(1,x,y,z);
        _v_out[index] = V.interior(1,x,y,z);
        _w_out[index] = W.interior(1,x,y,z);
        
        _txx_out[index] = Txx.interior(1,x,y,z);
        _tyy_out[index] = Tyy.interior(1,x,y,z);
        _tzz_out[index] = Tzz.interior(1,x,y,z);
         
    }
   }
  }
  
  {
    int dims[]={dimx, dimy, dimz};
    float spacing[]={h, h, h};
    opesci_dump_solution_vts("solution_pochoir", dims, spacing, _u_out,_v_out,_w_out,_txx_out,_tyy_out,_tzz_out);
  }
  
  {
    int dims[]={(int)round(sqrt(nrec)), (int)round(sqrt(nrec)), ntsteps};
    float spacing[]={h, h, dt};
    opesci_dump_receivers_vts("receivers_pochoir", dims, spacing, uss, vss, wss, pss);
  }

  return 0;
}
