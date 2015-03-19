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

#include <cassert>
#include <cstdlib>
#include <cmath>

#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

#include <vtkXMLStructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkZLibDataCompressor.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

// Abort execution.
void opesci_abort(std::string msg){
  std::cerr<<msg<<std::endl;
  exit(-1);
}

// Compute the field of Lame constants from subsurface model - i.e. the p-wave velocity, s-wave velocity and density fields.
// See http://scienceworld.wolfram.com/physics/LameConstants.html
void opesci_calculate_lame_costants(const std::vector<float> &vp, const std::vector<float> &vs, const std::vector<float> &rho,
				    std::vector<float> &mu, std::vector<float> &lam){
  size_t size=vp.size();
  assert(size==vs.size());
  assert(size==rho.size());
  mu.resize(size);
  lam.resize(size);
  
#pragma omp parallel for
  for(size_t i=0;i<size;i++){
    mu[i] = rho[i]*vs[i]*vs[i];
    lam[i] = rho[i]*(vp[i]*vp[i]-2.0*vs[i]*vs[i]);
  }
}

void opesci_dump_field_raw(std::string name, std::vector<float> &field){
  std::ofstream fh;
  fh.open (name+".raw", std::ios::out | std::ios::trunc | std::ios::binary); 
  fh.write((char *)field.data(), field.size()*sizeof(float));
  fh.close();
}

void opesci_dump_solution_vts(std::string name, const int dims[], const float spacing[],
			      std::vector<float> &u, std::vector<float> &v, std::vector<float> &w,
			      std::vector<float> &txx, std::vector<float> &tyy, std::vector<float> &tzz){

  vtkSmartPointer<vtkStructuredGrid> sg = vtkSmartPointer<vtkStructuredGrid>::New();
  sg->SetDimensions(dims[0], dims[1], dims[2]);

  {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(int k=0;k<dims[2];k++){
      float z = k*spacing[2];
      for(int j=0;j<dims[1];j++){
	float y = j*spacing[1];
	for(int i=0;i<dims[0];i++){
	  float x = i*spacing[0];
	  points->InsertNextPoint(x, y, z);
	}
      }
    }
    sg->SetPoints(points);
  }

  {
    vtkSmartPointer<vtkFloatArray> velocity = vtkSmartPointer<vtkFloatArray>::New();  
    velocity->SetName("velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);

    for(int k=0;k<dims[2];k++){
      for(int j=0;j<dims[1];j++){
	for(int i=0;i<dims[0];i++){
	  int index = k*dims[0]*dims[1]+j*dims[1]+i;
	  velocity->SetTuple3(index, u[index], v[index], w[index]);
	}
      }
    }
    sg->GetPointData()->AddArray(velocity);
  }

  {
    vtkSmartPointer<vtkFloatArray> pressure = vtkSmartPointer<vtkFloatArray>::New();  
    pressure->SetName("pressure");
    pressure->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
    pressure->SetNumberOfComponents(1);
    
    for(int k=0;k<dims[2];k++){
      for(int j=0;j<dims[1];j++){
	for(int i=0;i<dims[0];i++){
	  int index = k*dims[0]*dims[1]+j*dims[1]+i;
	  pressure->SetTuple1(index, (txx[index]+tyy[index]+tzz[index])/3);
	}
      }
    }
    sg->GetPointData()->AddArray(pressure);
  }

  vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(std::string(name+".vts").c_str());
  
  vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
  compressor->SetCompressionLevel(9);
  writer->SetCompressor(compressor);

  writer->SetInput(sg);
  writer->Write();
}

void opesci_dump_field_vts(std::string name, const int dims[], const float spacing[], std::vector<float> &field){
  assert(dims[0]*dims[1]*dims[2]==field.size());

  vtkSmartPointer<vtkStructuredGrid> sg = vtkSmartPointer<vtkStructuredGrid>::New();
  sg->SetDimensions(dims[0], dims[1], dims[2]);
  
  {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();  
    for(int k=0;k<dims[2];k++){
      float z = k*spacing[2];
      for(int j=0;j<dims[1];j++){
	float y = j*spacing[1];
	for(int i=0;i<dims[0];i++){
	  float x = i*spacing[0];
	  points->InsertNextPoint(x, y, z);
	}
      }
    }
    sg->SetPoints(points);
  }

  {
    vtkSmartPointer<vtkFloatArray> vtkfield = vtkSmartPointer<vtkFloatArray>::New();  
    vtkfield->SetName("field");
    vtkfield->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
    for(int k=0;k<dims[2];k++){
      for(int j=0;j<dims[1];j++){
	for(int i=0;i<dims[0];i++){
	  int index = k*dims[0]*dims[1]+j*dims[1]+i;
	  vtkfield->SetTuple1(index, field[index]);
	}
      }
    }
    sg->GetPointData()->AddArray(vtkfield);
  }

  vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(std::string(name+".vts").c_str());
  
  vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
  compressor->SetCompressionLevel(9);
  writer->SetCompressor(compressor);

  writer->SetInput(sg);
  writer->Write();
}

int opesci_read_simple_binary(const char *filename, std::vector<float> &array){
  std::ifstream infile(filename, std::ios::in | std::ios::binary);
  if(!infile.good()){
    std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open binary file "<<filename<<std::endl;
    return -1;
  }
  
  std::vector<unsigned char> buffer((std::istreambuf_iterator<char>(infile)),
				    std::istreambuf_iterator<char>());
    
  size_t size = buffer.size()/4;
  array.resize(size);
#pragma omp parallel for if (size >= 10000)
  for(size_t i=0;i<size;i++){
    array[i] = *((float*)&buffer[i*4]);
  }
  
  infile.close();
  
  return 0;
}

int opesci_get_souces(const char *xyz_filename, const char *xsrc_filename, const char *ysrc_filename, const char *zsrc_filename,
		      std::vector<float> &xyz_array, std::vector<float> &xsrc_array, std::vector<float> &ysrc_array, std::vector<float> &zsrc_array){
 std::ifstream infile(xyz_filename);
  if(!infile.good()){
    std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open source file "<<xyz_filename<<std::endl;
    return -1;
  }
  
  std::string line;
  std::getline(infile, line); // Read and ditch the header.
  while(!infile.eof()){
    std::getline(infile, line);
    if(line.empty())
      continue;
    
    std::stringstream ss(line);
    std::vector<float> xyz(3);
    ss>>xyz[0]>>xyz[1]>>xyz[2];
    xyz_array.insert(xyz_array.end(), xyz.begin(), xyz.end());
  }

  opesci_read_simple_binary(xsrc_filename, xsrc_array);
  opesci_read_simple_binary(ysrc_filename, ysrc_array);
  opesci_read_simple_binary(zsrc_filename, zsrc_array);
  
  return 0;
}

int opesci_get_receivers(const char *filename, std::vector<float> &array){
 std::ifstream infile(filename);
  if(!infile.good()){
    std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open receivers file "<<filename<<std::endl;
    return -1;
  }
  
  std::string line;
  std::getline(infile, line); // Read and ditch the header.
  while(!infile.eof()){
    std::getline(infile, line);
    if(line.empty())
      continue;
    
    std::stringstream ss(line);
    std::vector<float> xyz(3);
    ss>>xyz[0]>>xyz[1]>>xyz[2];
    array.insert(array.end(), xyz.begin(), xyz.end());
  }
  
  return 0;
}

float opesci_get_dt(const std::vector<float> &vp, float h){
  float maxv = 0;
  size_t size = vp.size();
#pragma omp parallel for reduction(max:maxv)
  for(size_t i=0;i<size;++i){
    if(vp[i]>maxv){
      maxv = vp[i];
    }
  }

  return (6.0/7.0)*(1./sqrt(3.0))*(h/maxv);
}

/* 
 * Discrete Fourier transform
 * by Project Nayuki, 2014. Public domain.
 * http://www.nayuki.io/page/how-to-implement-the-discrete-fourier-transform
 * 
 * Computes the discrete Fourier transform (DFT) of the given vector.
 * All the array arguments must have the same length.
 */

void compute_dft(const float inreal[], const float inimag[], float outreal[], float outimag[], int n){
#pragma omp parallel for
  for(int k=0;k<n;++k) {  /* For each output element */
    float sumreal = 0;
    float sumimag = 0;
    for(int t=0;t<n;++t) {  /* For each input element */
      float angle = 2 * M_PI * t * k / n;
      sumreal +=  inreal[t] * cos(angle) + inimag[t] * sin(angle);
      sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
    }
    outreal[k] = sumreal;
    outimag[k] = sumimag;
  }
}

/* 
 * If dt < sdt, the source is resampled by perfoming a discrete
 * fourier transfor, filling the resulting spectrum with as many
 * zeroes as needed and inverse fourier transforming to get the
 * desired sampling without changing the frequency content.
 *
 * If sdt < dt, the source is resampled as above but removing samples
 * of the spectrum until the desired sample rate is desired.
 */
void resample(const std::vector<float> &src, float dt, double sdt, std::vector<float> &resampled_src){
  // Break out early if necessary.
  if(fabs(dt-sdt)<std::numeric_limits<float>::epsilon()*(dt+sdt)){
    resampled_src = src;
    return;
  }
  
  // Calculate new source size.
  int snt = src.size();
  int snt2 = (int)round(snt*sdt/dt);

  std::vector<float> imsrc(snt, 0.0);
  
  std::vector<float> ft_resrc(snt);
  std::vector<float> ft_imsrc(snt);

  compute_dft(src.data(), imsrc.data(), ft_resrc.data(), ft_imsrc.data(), snt);
  
  // Normalize
  float nrm = 1./sqrt(snt);
  for(auto &i : ft_resrc)
    i*=nrm;
  for(auto &i : ft_imsrc)
    i*=nrm;

  std::vector<float> ft_resrc2(snt2, 0.0);
  std::vector<float> ft_imsrc2(snt2, 0.0);
  
  if(dt<sdt){
    // Add zeroes in the middle of the real and imaginary part of the
    // transform. This is the same as adding zeroes at the end of the
    // spectrum.
    int midpoint = snt/2;
    for(int i=0;i<midpoint;i++){
      ft_resrc2[i] = ft_resrc[i];
      ft_imsrc2[i] = ft_imsrc[i];
    }
    int offset = (snt2-snt);
    for(int i=midpoint;i<snt;i++){
      ft_resrc2[offset+i] = ft_resrc[i];
      ft_imsrc2[offset+i] = ft_imsrc[i];
    }
  }else if(dt>sdt){
    // Substract the middle of the real and imaginary part of the transform.
    int midpoint = snt2/2;
    for(int i=0;i<midpoint;i++){
      ft_resrc2[i] = ft_resrc[i];
      ft_imsrc2[i] = ft_imsrc[i];
    }
    int offset = (snt-snt2);
    for(int i=midpoint;i<snt2;i++){
      ft_resrc2[i] = ft_resrc[offset+i];
      ft_imsrc2[i] = ft_imsrc[offset+i];
    }
  }

  // Inverse fourier transform.
  for(auto &i : ft_imsrc2)
    i*=-1;
  
  resampled_src.resize(snt2);
  std::vector<float> imsrc2(snt2);
  
  compute_dft(ft_resrc2.data(), ft_imsrc2.data(), resampled_src.data(), imsrc2.data(), snt2);
  
  // Normalise
  float nrm2 = 1.0/sqrt((float)snt2);
  for(auto &i : resampled_src)
    i*=nrm2;
}
// #include <fenv.h> 

int main(){
  // feenableexcept(FE_INVALID | FE_OVERFLOW);

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
  if(opesci_get_souces(geomsrc.c_str(), xsrcname.c_str(), ysrcname.c_str(), zsrcname.c_str(), coorsrc, xsrc, ysrc, zsrc)){
    opesci_abort("Cannot read sources.\n");
  }
  int nsrc = coorsrc.size()/3;

  // Get receivers.
  std::vector<float> coorrec;
  if(opesci_get_receivers(geomrec.c_str(), coorrec)){
    opesci_abort("Cannot read receivers.\n");
  }
  int nrec = coorrec.size()/3;

  float dt = opesci_get_dt(vp, h);
  int ntsteps = (int)(maxt/dt);
  
  // Resample sources if required.
  assert(nsrc==1);
  std::vector<float> resampled_src;
  resample(xsrc, dt, sdt, resampled_src);
  xsrc.swap(resampled_src);

  resample(ysrc, dt, sdt, resampled_src);
  ysrc.swap(resampled_src);

  resample(zsrc, dt, sdt, resampled_src);
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
