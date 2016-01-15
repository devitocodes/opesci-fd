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
#include <cstring>

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

void opesci_dump_field_raw(std::string name, std::vector<float> &field)
{
    std::ofstream fh;
    fh.open (name+".raw", std::ios::out | std::ios::trunc | std::ios::binary);
    fh.write((char *)field.data(), field.size()*sizeof(float));
    fh.close();
}

void opesci_dump_solution_vts(std::string name, const int dims[], const float spacing[],
                              std::vector<float> &u, std::vector<float> &v, std::vector<float> &w,
                              std::vector<float> &txx, std::vector<float> &tyy, std::vector<float> &tzz)
{
#ifdef VTK_FOUND
    vtkSmartPointer<vtkStructuredGrid> sg = vtkSmartPointer<vtkStructuredGrid>::New();
    sg->SetDimensions(dims[0], dims[1], dims[2]);

    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for(int k=0; k<dims[2]; k++) {
            float z = k*spacing[2];
            for(int j=0; j<dims[1]; j++) {
                float y = j*spacing[1];
                for(int i=0; i<dims[0]; i++) {
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

        for(int k=0; k<dims[2]; k++) {
            for(int j=0; j<dims[1]; j++) {
                for(int i=0; i<dims[0]; i++) {
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

        for(int k=0; k<dims[2]; k++) {
            for(int j=0; j<dims[1]; j++) {
                for(int i=0; i<dims[0]; i++) {
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

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(sg);
#else
    writer->SetInputData(sg);
#endif
    writer->Write();
#else
    opesci_abort("ERROR: OPESCI built without VTK support. Cannot dump VTK files.");
#endif
}

void opesci_dump_field_vts(std::string name, const int dims[], const float spacing[], std::vector<float> &field)
{
#ifdef VTK_FOUND
    assert(dims[0]*dims[1]*dims[2]==field.size());

    vtkSmartPointer<vtkStructuredGrid> sg = vtkSmartPointer<vtkStructuredGrid>::New();
    sg->SetDimensions(dims[0], dims[1], dims[2]);

    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(dims[0]*dims[1]*dims[2]);

        vtkSmartPointer<vtkFloatArray> vtkfield = vtkSmartPointer<vtkFloatArray>::New();
        vtkfield->SetName("field");
        vtkfield->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);

        vtkIdType pcnt=0;
        for(int k=0; k<dims[2]; k++) {
            float z = k*spacing[2];
            for(int j=0; j<dims[1]; j++) {
                float y = j*spacing[1];
                for(int i=0; i<dims[0]; i++) {
                    float x = i*spacing[0];

                    points->SetPoint(pcnt, x, y, z);
                    vtkfield->SetTuple1(pcnt, field[pcnt]);

                    pcnt++;
                }
            }
        }
        sg->SetPoints(points);
        sg->GetPointData()->AddArray(vtkfield);
    }

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(std::string(name+".vts").c_str());

    vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
    compressor->SetCompressionLevel(1);
    writer->SetCompressor(compressor);

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(sg);
#else
    writer->SetInputData(sg);
#endif
    writer->Write();
#else
    opesci_abort("ERROR: OPESCI built without VTK support. Cannot dump VTK files.");
#endif
}

void opesci_dump_receivers_vts(std::string name, const int dims[], const float spacing[],
                               std::vector<float> &uss, std::vector<float> &vss, std::vector<float> &wss, std::vector<float> &pss)
{
#ifdef VTK_FOUND
    vtkSmartPointer<vtkStructuredGrid> sg = vtkSmartPointer<vtkStructuredGrid>::New();
    sg->SetDimensions(dims[0], dims[1], dims[2]);

    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for(int k=0; k<dims[2]; k++) {
            float z = k*spacing[2];
            for(int j=0; j<dims[1]; j++) {
                float y = j*spacing[1];
                for(int i=0; i<dims[0]; i++) {
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
        for(int k=0; k<dims[2]; k++) {
            for(int j=0; j<dims[1]; j++) {
                for(int i=0; i<dims[0]; i++) {
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
        for(int k=0; k<dims[2]; k++) {
            for(int j=0; j<dims[1]; j++) {
                for(int i=0; i<dims[0]; i++) {
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
        for(int k=0; k<dims[2]; k++) {
            for(int j=0; j<dims[1]; j++) {
                for(int i=0; i<dims[0]; i++) {
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
        for(int k=0; k<dims[2]; k++) {
            for(int j=0; j<dims[1]; j++) {
                for(int i=0; i<dims[0]; i++) {
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

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(sg);
#else
    writer->SetInputData(sg);
#endif
    writer->Write();
#else
    opesci_abort("ERROR: OPESCI built without VTK support. Cannot dump VTK files.");
#endif
}

int opesci_read_simple_binary(const char *filename, std::vector<float> &array)
{
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if(!infile.good()) {
        std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open binary file "<<filename<<std::endl;
        return -1;
    }

    std::vector<unsigned char> buffer((std::istreambuf_iterator<char>(infile)),
                                      std::istreambuf_iterator<char>());

    size_t size = buffer.size()/4;
    array.resize(size);
    #pragma omp parallel for if (size >= 10000)
    for(size_t i=0; i<size; i++) {
        array[i] = *((float*)&buffer[i*4]);
    }

    infile.close();

    return 0;
}

int opesci_read_simple_binary_ptr(const char *filename, float *array, int size)
{
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if(!infile.good()) {
        std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open binary file "<<filename<<std::endl;
        return -1;
    }

    std::vector<unsigned char> buffer((std::istreambuf_iterator<char>(infile)),
                                      std::istreambuf_iterator<char>());

    size_t filesize = buffer.size()/4;
    if (filesize>size) {
        std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Input file "<<filename<<" size larger than array size "<<std::endl;
    }
    #pragma omp parallel for if (size >= 10000)
    for(size_t i=0; i<size; i++) {
        array[i] = *((float*)&buffer[i*4]);
    }

    infile.close();

    return 0;
}


int opesci_read_souces(const char *xyz_filename, const char *xsrc_filename, const char *ysrc_filename, const char *zsrc_filename,
                       std::vector<float> &xyz_array, std::vector<float> &xsrc_array, std::vector<float> &ysrc_array, std::vector<float> &zsrc_array)
{
    std::ifstream infile(xyz_filename);
    if(!infile.good()) {
        std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open source file "<<xyz_filename<<std::endl;
        return -1;
    }

    std::string line;
    std::getline(infile, line); // Read and ditch the header.
    while(!infile.eof()) {
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

int opesci_read_receivers(const char *filename, std::vector<float> &array)
{
    std::ifstream infile(filename);
    if(!infile.good()) {
        std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open receivers file "<<filename<<std::endl;
        return -1;
    }

    std::string line;
    std::getline(infile, line); // Read and ditch the header.
    while(!infile.eof()) {
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

// Lifted from segy_util.c
/* convert IBM FORTRAN REAL*4 to IEEE float */
float real2float(const char *xreal)
{
    union {
        unsigned long l;
        char c[4];
    } yreal;
    long vorzeichen, exponent, mantisse, i;

    for (i=0; i<4; i++)
        yreal.c[i]=xreal[i];
    vorzeichen=((yreal.l & 0x80000000) != 0);
    exponent=((yreal.l & 0x7f000000) >> 24)-64;
    mantisse=(yreal.l & 0x00ffffff);
    if (vorzeichen)
        return(-mantisse/16777216.0*pow(16.0,(double)exponent));
    else
        return( mantisse/16777216.0*pow(16.0,(double)exponent));
}

// Unpack required data from header. Start by assuming endian swap is not required.
inline int16_t unpack_int16(const char *target, bool swap_endian)
{
    if(swap_endian) {
        char buffer[2] = {target[1], target[0]};
        return *(int16_t *)(buffer);
    } else {
        return *(int16_t *)(target);
    }
}

inline int32_t unpack_int32(const char *target, bool swap_endian)
{
    if(swap_endian) {
        char buffer[4] = {target[3], target[2], target[1], target[0]};
        return *(int32_t *)(buffer);
    } else {
        return *(int32_t *)(target);
    }
}

inline float unpack_ibmfloat(const char *target, bool swap_endian)
{
    if(swap_endian) {
        char buffer[4] = {target[3], target[2], target[1], target[0]};
        return real2float(buffer);
    } else {
        return real2float(target);
    }
}

// Assuming http://www.seg.org/documents/10161/77915/seg_y_rev1.pdf
int opesci_read_model_segy(const char *filename, std::vector<float> &array, int dim[], float spacing[])
{
    std::ifstream infile(filename);
    if(!infile.good()) {
        std::cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Failed to open SEG-Y file "<<filename<<std::endl;
        return -1;
    }

    // Get the size of the file.
    infile.seekg (0, infile.end);
    long filesize = infile.tellg();

    // Read in header
    char header_buffer[3600];
    infile.seekg(0, infile.beg);
    infile.read(header_buffer, 3600);

    bool swap_endian=false;
    int ntraces, tracesize, format_code;
    int Nx, Ny, Nz;
    for(int i=0; i<2; i++) { // Enables a re-try of header read.
        // 3201 - 3204 	Job identification number.
        // std::cout<<"3201 - 3204 	Job identification number: "<<unpack_int32(header_buffer+3200, swap_endian)<<std::endl;

        // 3205 - 3208 	Line number.
        // std::cout<<"3205 - 3208 	Line number: "<<unpack_int32(header_buffer+3204, swap_endian)<<std::endl;

        // 3209 - 3212 	Reel number.
        // std::cout<<"3209 - 3212 	Reel number: "<<unpack_int32(header_buffer+3208, swap_endian)<<std::endl;

        // 3213 - 3214 	Number of data traces per record.
        Nx = unpack_int16(header_buffer+3212, swap_endian);

        // 3215 - 3216 	Number of auxiliary traces per record.
        // std::cout<<"3215 - 3216 	Number of auxiliary traces per record: "<<unpack_int16(header_buffer+3214, swap_endian)<<std::endl;

        // 3217 - 3218 	Sample interval, microseconds, this file (reel).
        // std::cout<<"3217 - 3218 	Sample interval, microseconds, this file (reel): "<<unpack_int16(header_buffer+3216, swap_endian)<<std::endl;

        // 3219 - 3220 	Sample interval, microseconds, original field recording.
        // 3221 - 3222 	Number of samples per data trace, this file (reel).
        Nz = unpack_int16(header_buffer+3220, swap_endian);

        // 3223 - 3224 	Number of samples per data trace, original field recording.
        // 3225 - 3226 	Data sample format code: 1 = 4-byte IBM floating-point
        //                                           2 = 4-byte, two's complement integer
        //                                           3 = 2-byte, two's complement integer
        //                                           4 = 4-byte fixed-point with gain (obsolete)
        //                                           5 = 4-byte IEEE floating-point
        //                                           6 = Not currently used
        //                                           7 = Not currently used
        //                                           8 = 1-byte, two's complement integer
        format_code = unpack_int16(header_buffer+3224, swap_endian);
        if(format_code<1 || format_code>8) {
            swap_endian = true;
            format_code = unpack_int16(header_buffer+3224, swap_endian);

            // Try again...
            if(format_code<1 || format_code>8) {
                std::cerr<<"ERROR: unsupported data sample format code "<<format_code<<std::endl;
                return -1;
            }

            continue; // restart header read
        }

        // 3227 - 3228 	CDP fold.
        // std::cout<<"3227 - 3228 	CDP fold: "<<unpack_int16(header_buffer+3226, swap_endian)<<std::endl;

        // 3229 - 3230 	Trace sorting code: Trace sorting code (i.e. type of ensemble) :
        //                                      -1 = Other (should be explained in user Extended Textual File Header stanza
        //                                       0 = Unknown
        //                                       1 = As recorded (no sorting)
        //                                       2 = CDP ensemble
        //                                       3 = Single fold continuous profile
        //                                       4 = Horizontally stacked
        //                                       5 = Common source point
        //                                       6 = Common receiver point
        //                                       7 = Common offset point
        //                                       8 = Common mid-point
        //                                       9 = Common conversion point
        // std::cout<<"3229 - 3230 	Trace sorting code: "<<unpack_int16(header_buffer+2328, swap_endian)<<std::endl;

        // 3231 - 3232 	Vertical sum code: 1 = no sum 2 = two sum ... N = N sum (N = 32,767)
        // std::cout<<"3231 - 3232 	Vertical sum code: "<<unpack_int16(header_buffer+2330, swap_endian)<<std::endl;

        // 3233 - 3234 	Sweep frequency at start.
        // 3235 - 3236 	Sweep frequency at end.
        // 3237 - 3238 	Sweep length, ms.
        // 3239 - 3240 	Sweep type code: 1 = linear
        //                                   2 = parabolic
        //                                   3 = exponential
        //                                   4 = other
        // 3241 - 3242 	Trace number of sweep channel.
        // std::cout<<"3241 - 3242 	Trace number of sweep channel: "<<unpack_int16(header_buffer+2340, swap_endian)<<std::endl;

        // 3243 - 3244 	Sweep trace taper length, ms, at start if tapered.
        // std::cout<<"3243 - 3244 	Sweep trace taper length, ms, at start if tapered: "<<unpack_int16(header_buffer+2342, swap_endian)<<std::endl;

        // 3245 - 3246 	Sweep trace taper length, ms, at end.
        // std::cout<<"3245 - 3246 	Sweep trace taper length, ms, at end: "<<unpack_int16(header_buffer+2344, swap_endian)<<std::endl;

        // 3247 - 3248 	Taper type: 1 = linear 2 = cos 3 = other
        // std::cout<<"3247 - 3248 	Taper type: 1 = linear 2 = cos 3 = other: "<<unpack_int16(header_buffer+2346, swap_endian)<<std::endl;

        // 3249 - 3250 	Correlated data traces: 1 = no 2 = yes
        // std::cout<<"3249 - 3250 	Correlated data traces: 1 = no 2 = yes"<<unpack_int16(header_buffer+2348, swap_endian)<<std::endl;

        // 3251 - 3252 	Binary gain recovered: 1 = yes 2 = no
        // std::cout<<"3251 - 3252 	Binary gain recovered: 1 = yes 2 = no: "<<unpack_int16(header_buffer+2350, swap_endian)<<std::endl;

        // 3253 - 3254 	Amplitude recovery method: 1 = none 2 = spherical divergence 3 = AGC 4 = other
        // std::cout<<"3253 - 3254 	Amplitude recovery method: 1 = none 2 = spherical divergence 3 = AGC 4 = other: "<<unpack_int16(header_buffer+2352, swap_endian)<<std::endl;

        // 3255 - 3256 	Measurement system: 1 = meters 2 = feet
        // std::cout<<"3255 - 3256 	Measurement system: 1 = meters 2 = feet: "<<unpack_int16(header_buffer+2354, swap_endian)<<std::endl;

        // 3257 - 3258 	Impulse signal: 1 = Upward = negative number.
        //                                  2 = Upward = positive number.
        // std::cout<<"3257 - 3258 	Impulse signal: "<<unpack_int16(header_buffer+2356, swap_endian)<<std::endl;

        // 3259 - 3260 	Vibratory polarity code - seismic signal lags pilot signal by: 1 = 337.5 - 22.5 degrees
        //                                                                                 2 = 22.5 - 67.5 degrees
        //                                                                                 3 = 67.5 - 112.5 degrees
        //                                                                                 4 = 112.5 - 157.5 degrees
        //                                                                                 5 = 157.5 - 202.5 degrees
        //                                                                                 6 = 202.5 - 247.5 degrees
        // std::cout<<"3259 - 3260 	Vibratory polarity code - seismic signal lags pilot signal by: "<<unpack_int16(header_buffer+2358, swap_endian)<<std::endl;

        // Set tracesize
        tracesize = 240+4*Nz;

        // Get number of traces.
        ntraces = (filesize-3600)/tracesize;

        Ny = ntraces/Nx;

        break;
    }

    dim[0] = Nx;
    dim[1] = Ny;
    dim[2] = Nz;

    int xysize=dim[0]*dim[1];
    int zysize=dim[2]*dim[1];

    array.resize(dim[0]*dim[1]*dim[2]);
    if(format_code==1) {
        char trace_buffer[tracesize];
        float x0[2], x1[2], scale_xy;
        for(int i=0; i<ntraces; i++) {
            // See Trace header in http://www.seg.org/documents/10161/77915/seg_y_rev1.pdf
            infile.read(trace_buffer, tracesize);

            int ix = i%Nx;
            int iy = i/Nx;

            if(i==0) {
                scale_xy = unpack_int16(trace_buffer+70, swap_endian);
                if(scale_xy<0)
                    scale_xy = 1.0/fabs(scale_xy);
                x0[0] = scale_xy*unpack_int32(trace_buffer+72, swap_endian);
                x0[1] = scale_xy*unpack_int32(trace_buffer+76, swap_endian);
            } else if(i==1) {
                x1[0] = scale_xy*unpack_int32(trace_buffer+72, swap_endian);
                x1[1] = scale_xy*unpack_int32(trace_buffer+76, swap_endian);

                float dx = std::max(fabs(x0[0]-x1[0]), fabs(x0[1]-x1[1]));

                spacing[0] = dx;
                spacing[1] = dx;
                spacing[2] = scale_xy*unpack_int16(trace_buffer+116, swap_endian); // For model data this is overloaded as meters.
            }
            assert(unpack_int16(trace_buffer+114, swap_endian)==Nz);

#ifndef NDEBUG
            if(i>0) {
                int _ix = (scale_xy*unpack_int32(trace_buffer+72, swap_endian)-x0[0])/spacing[0] + 0.5; // +0.5 to protext against round-off
                int _iy = (scale_xy*unpack_int32(trace_buffer+76, swap_endian)-x0[1])/spacing[1] + 0.5;
                assert(ix==_ix);
                assert(iy==_iy);
            }
#endif

            for(int iz=0; iz<Nz; iz++) {
                char *next = trace_buffer+240+iz*4;
                array[ix+iy*Nx+iz*xysize] = unpack_ibmfloat(next, swap_endian);
            }
        }
    } else {
        std::cerr<<"ERROR: format code "<<format_code<<" not yet supported"<<std::endl;
        return -1;
    }

    return 0;
}

void opesci_dump_field_vts_3d(std::string name, const int dims[], const float spacing[], int margin, float *field)
{
#ifdef VTK_FOUND

    vtkSmartPointer<vtkStructuredGrid> sg = vtkSmartPointer<vtkStructuredGrid>::New();
    sg->SetDimensions(dims[0], dims[1], dims[2]);

    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for(int i=0; i<dims[0]; i++) {
            float x = (i-margin)*spacing[0];
            for(int j=0; j<dims[1]; j++) {
                float y = (j-margin)*spacing[1];
                for(int k=0; k<dims[2]; k++) {
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
        for(int i=0; i<dims[0]; i++) {
            for(int j=0; j<dims[1]; j++) {
                for(int k=0; k<dims[2]; k++) {
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

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(sg);
#else
    writer->SetInputData(sg);
#endif
    writer->Write();
#else
    opesci_abort("ERROR: OPESCI built without VTK support. Cannot dump VTK files.");
#endif
}
