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

  int index = 0;
  
  
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
  Pochoir_Array_3D(float, 1) U(dimx, dimy, dimz), V(dimx, dimy, dimz), W(dimx, dimy, dimz),
    Txx(dimx, dimy, dimz), Tyy(dimx, dimy, dimz), Tzz(dimx, dimy, dimz),
    Tyz(dimx, dimy, dimz), Txz(dimx, dimy, dimz), Txy(dimx, dimy, dimz);

  // Subsurface model
  Pochoir_Array_3D(float, 1) beta(dimx, dimy, dimz), lambda(dimx, dimy, dimz), mu(dimx, dimy, dimz);
  
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
  fd_3D.Register_Array(U);
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
