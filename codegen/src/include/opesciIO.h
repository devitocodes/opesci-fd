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

#ifndef OPESCIIO_H
#define OPESCIIO_H
#include <string>
#include <vector>

void opesci_dump_field_raw(std::string name, std::vector<float> &field);
void opesci_dump_solution_vts(std::string name, const int dims[], const float spacing[],
			      std::vector<float> &u, std::vector<float> &v, std::vector<float> &w,
			      std::vector<float> &txx, std::vector<float> &tyy, std::vector<float> &tzz);
void opesci_dump_field_vts(std::string name, const int dims[], const float spacing[], std::vector<float> &field);
void opesci_dump_receivers_vts(std::string name, const int dims[], const float spacing[],
			       std::vector<float> &uss, std::vector<float> &vss, std::vector<float> &wss, std::vector<float> &pss);

int opesci_read_simple_binary(const char *filename, std::vector<float> &array);
int opesci_read_simple_binary_ptr(const char *filename, float *array);
int opesci_read_souces(const char *xyz_filename, const char *xsrc_filename, const char *ysrc_filename, const char *zsrc_filename,
		       std::vector<float> &xyz_array, std::vector<float> &xsrc_array, std::vector<float> &ysrc_array, std::vector<float> &zsrc_array);
int opesci_read_receivers(const char *filename, std::vector<float> &array);
int opesci_read_model_segy(const char *filename, std::vector<float> &array, int dim[], float spacing[]);

void opesci_dump_field_vts_3d(std::string name, const int dims[], const float spacing[], int margin, float *field);

#endif
