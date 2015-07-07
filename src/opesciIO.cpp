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

#include <algorithm>
#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

#ifdef VTK_FOUND
#include <vtkXMLStructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkZLibDataCompressor.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#endif

#include "opesciIO.h"
#include "opesciHandy.h"

void opesci_dump_field_raw(std::string name, std::vector<float> &field){
  std::ofstream fh;
  fh.open (name+".raw", std::ios::out | std::ios::trunc | std::ios::binary); 
  fh.write((char *)field.data(), field.size()*sizeof(float));
  fh.close();
}

void opesci_dump_solution_vts(std::string name, const int dims[], const float spacing[],
			      std::vector<float> &u, std::vector<float> &v, std::vector<float> &w,
			      std::vector<float> &txx, std::vector<float> &tyy, std::vector<float> &tzz){
#ifdef VTK_FOUND
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
#else
  opesci_abort("ERROR: OPESCI built without VTK support. Cannot dump VTK files.");
#endif
}

void opesci_dump_field_vts(std::string name, const int dims[], const float spacing[], std::vector<float> &field){
#ifdef VTK_FOUND
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
	  int index = k*dims[0]*dims[1]+j*dims[0]+i;
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
#else
  opesci_abort("ERROR: OPESCI built without VTK support. Cannot dump VTK files.");
#endif
}

void opesci_dump_receivers_vts(std::string name, const int dims[], const float spacing[],
			       std::vector<float> &uss, std::vector<float> &vss, std::vector<float> &wss, std::vector<float> &pss){
#ifdef VTK_FOUND
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
    vtkfield->SetName("vss");
    vtkfield->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
    for(int k=0;k<dims[2];k++){
      for(int j=0;j<dims[1];j++){
	for(int i=0;i<dims[0];i++){
	  int index = k*dims[0]*dims[1]+j*dims[1]+i;
	  vtkfield->SetTuple1(index, vss[index]);
	}
      }
    }
    sg->GetPointData()->AddArray(vtkfield);
  }

  {
    vtkSmartPointer<vtkFloatArray> vtkfield = vtkSmartPointer<vtkFloatArray>::New();  
    vtkfield->SetName("wss");
    vtkfield->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
    for(int k=0;k<dims[2];k++){
      for(int j=0;j<dims[1];j++){
	for(int i=0;i<dims[0];i++){
	  int index = k*dims[0]*dims[1]+j*dims[1]+i;
	  vtkfield->SetTuple1(index, wss[index]);
	}
      }
    }
    sg->GetPointData()->AddArray(vtkfield);
  }

  {
    vtkSmartPointer<vtkFloatArray> vtkfield = vtkSmartPointer<vtkFloatArray>::New();  
    vtkfield->SetName("pss");
    vtkfield->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
    for(int k=0;k<dims[2];k++){
      for(int j=0;j<dims[1];j++){
	for(int i=0;i<dims[0];i++){
	  int index = k*dims[0]*dims[1]+j*dims[1]+i;
	  vtkfield->SetTuple1(index, pss[index]);
	}
      }
    }
    sg->GetPointData()->AddArray(vtkfield);
  }

  {
    vtkSmartPointer<vtkFloatArray> vtkfield = vtkSmartPointer<vtkFloatArray>::New();  
    vtkfield->SetName("uss");
    vtkfield->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
    for(int k=0;k<dims[2];k++){
      for(int j=0;j<dims[1];j++){
	for(int i=0;i<dims[0];i++){
	  int index = k*dims[0]*dims[1]+j*dims[1]+i;
	  vtkfield->SetTuple1(index, uss[index]);
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
#else
  opesci_abort("ERROR: OPESCI built without VTK support. Cannot dump VTK files.");
#endif
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

int opesci_read_souces(const char *xyz_filename, const char *xsrc_filename, const char *ysrc_filename, const char *zsrc_filename,
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

int opesci_read_receivers(const char *filename, std::vector<float> &array){
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

int opesci_read_model_segy(const char *filename, std::vector<float> &array, int dim[], float spacing[]){
  std::ifstream infile(filename);
  if(!infile.good()){
    std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open SEG-Y file "<<filename<<std::endl;
    return -1;
  }

  // Get the size of the file.
  infile.seekg (0, infile.end);
  long filesize = infile.tellg();
 
  // Get the number of depth layers.
  infile.seekg(3220, std::ios::beg);
  char ns[2];
  infile.read(ns, 2);
  dim[0] = *(int16_t *)(ns);
  
  // Calculate number of traces in this file and the other dimensions
  // of the dataset. Include a heuristic for checking endian.
  bool swap_endian=false;
  int tracesize = 240+4*dim[0];
  int ntraces = (filesize-3600)/tracesize + 0.5;
  dim[1] = (int)(sqrt(ntraces)+0.5);
  dim[2] = dim[1];
  if(dim[1]*dim[2]!=ntraces){
    swap_endian=true;
    std::swap(ns[0], ns[1]);
    dim[0] = *(int16_t *)(ns);
    
    tracesize = 240+4*dim[0];
    ntraces = (filesize-3600)/tracesize + 0.5;
    dim[1] = (int)(sqrt(ntraces)+0.5);
    dim[2] = dim[1];
  }

  array.resize(dim[0]*dim[1]*dim[2]);  
  if(swap_endian){
    char buffer[dim[0]*4];
    for(int i=0;i<ntraces;i++){
      infile.seekg(3600+i*tracesize+240, std::ios::beg);
      infile.read(buffer, dim[0]*4);
      
      for(int j=0;j<dim[0];j++){
	std::swap(buffer[j*4], buffer[j*4+3]);
	std::swap(buffer[j*4+1], buffer[j*4+2]);
      }
      memcpy(array.data()+i*dim[0], buffer, dim[0]*4);
    }
  }else{
    for(int i=0;i<ntraces;i++){
      infile.seekg(3600+i*tracesize+240, std::ios::beg);
      infile.read((char *)(array.data()+i*dim[0]), dim[0]*4);
    }
  }

  return 0;
}



void opesci_dump_field_vts_3d(std::string name, const int dims[], const float spacing[], int margin, float *field){
#ifdef VTK_FOUND

  vtkSmartPointer<vtkStructuredGrid> sg = vtkSmartPointer<vtkStructuredGrid>::New();
  sg->SetDimensions(dims[0], dims[1], dims[2]);
  
  {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(int i=0;i<dims[0];i++){
      float x = (i-margin)*spacing[0];
      for(int j=0;j<dims[1];j++){
        float y = (j-margin)*spacing[1];
        for(int k=0;k<dims[2];k++){
          float z = (k-margin)*spacing[2];
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
    for(int i=0;i<dims[0];i++){
      for(int j=0;j<dims[1];j++){
        for(int k=0;k<dims[2];k++){
          int index = i*dims[1]*dims[2]+j*dims[2]+k;
          vtkfield->SetTuple1(index, *(field+index));
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
#else
  opesci_abort("ERROR: OPESCI built without VTK support. Cannot dump VTK files.");
#endif
}