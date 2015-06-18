<%include file="copyright.txt"/>

<%include file="common_include.txt"/>

#include <pochoir.hpp>
#include <math.h>

<%include file="pochoir_bc_2d.txt"/>

float h = 0.25;
int dimx = int dimy = (int) 1/h;
float lambda = 0.5, mu = 0.25; // lame constatns
float rho = 1.0; // density
float _tmax = 5.0; // simulate until
float dt = 0.5*h /(2*Vp)
int _ntsteps = (int) _tmax/dt;
float Vp = sqrt((lambda + 2*mu)/rho)

int main(){

  Pochoir_Shape_2D fd_shape_2D[] = {
    {1,0,0},
    {0,1,0},
    {0,0,1},
    {0,2,0},
    {0,0,2},
    {0,-1,0},
    {0,0,-1},
    {0,-2,0},
    {0,0,-2}};
  
  // Prognostic fields.
  Pochoir_Array<float, 3> U(dimx, dimy), V(dimx, dimy), Txx(dimx, dimy), Tyy(dimx, dimy), Tyz(dimx, dimy);

  Pochoir_2D fd_2D(fd_shape_2D);
  
  U.Register_Boundary(fd_bv_2D_U);
  V.Register_Boundary(fd_bv_2D_V);
  Txx.Register_Boundary(fd_bv_2D_Txx);
  Tyy.Register_Boundary(fd_bv_2D_Tyy);
  Txy.Register_Boundary(fd_bv_2D_Txy);

  fd_2D.Register_Array(U);
  fd_2D.Register_Array(V);
  fd_2D.Register_Array(Txx);
  fd_2D.Register_Array(Tyy);
  fd_2D.Register_Array(Txy);

  //Pochoir_Domain X(0, dimx-1), Y(0, dimy-1);
  //fd_2D.Register_Domain(X, Y);
  //U.Register_Shape(fd_shape_3D);
  //V.Register_Shape(fd_shape_3D);
  
  // Initialization of prognostic fields
  for(int i=0;i<dimx;i++){
    for(int j=0;j<dimy;j++){
      float x = i*h;
      float y = j*h;
      U(0,x,y) = ${U_init};
      V(0,x,y) = ${V_init};

      Txx(0,x,y) = ${Txx_init};
      Tyy(0,x,y) = ${Tyy_init};
      Txy(0,x,y) = ${Txy_init};
    }
  }
  
  
  Pochoir_Kernel_3D(fd_3D_velocity, t, x, y, z)
    // Update velocity
    ${U};
    ${V};
    ${W};
           
  Pochoir_Kernel_End

  Pochoir_Kernel_3D(fd_3D_stress, t, x, y, z)
    // Update stress
    ${Txx};
    ${Tyy};
    ${Tzz};
    ${Tyz};
    ${Txz};
    ${Txy};

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
            Txx(0,_sx,_sy,_sz) = Txx(1,_sx,_sy,_sz) -= xsrc[_t/2]/3;
            Tyy(0,_sx,_sy,_sz) = Tyy(1,_sx,_sy,_sz) -= ysrc[_t/2]/3;
            Tzz(0,_sx,_sy,_sz) = Tzz(1,_sx,_sy,_sz) -= zsrc[_t/2]/3;
          }
        }
      } //end of parallel region
    } //end else
  } //end for
  
  for(int x=0;x<dimx;++x){
    for(int y=0;y<dimy;++y){
      for(int z=0;z<dimz;++z){
       
        int index = x*dimy*dimz+y*dimz+z;
        
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
    opesci_dump_solution_vts("solution_pochoir_tmpl", dims, spacing, _u_out,_v_out,_w_out,_txx_out,_tyy_out,_tzz_out);
  }
  
  {
    int dims[]={(int)round(sqrt(nrec)), (int)round(sqrt(nrec)), ntsteps};
    float spacing[]={h, h, dt};
    opesci_dump_receivers_vts("receivers_pochoir_tmpl", dims, spacing, uss, vss, wss, pss);
  }

  return 0;
}
