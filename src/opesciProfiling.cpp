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
#include "opesciProfiling.h"
#if defined(OPESCI_HAVE_PAPI)
#include "papi.h"
#endif
#include "omp.h"

#define OPESCI_PAPI_WARN "WARNING:: PAPI error: Flops/s counters are not reliable!\n"
#define OPESCI_PAPI_MISSING "WARNING:: PAPI not found: Please re-build Opesci-FD with PAPI libraries\n"

int opesci_papi_init()
{
#if defined(OPESCI_HAVE_PAPI)
    int err, version;
    version = PAPI_library_init(PAPI_VER_CURRENT);
    if (version != PAPI_VER_CURRENT) printf(OPESCI_PAPI_WARN);

    err = PAPI_thread_init((unsigned long (*)()) omp_get_thread_num);
    if (err != PAPI_OK) printf(OPESCI_PAPI_WARN);
    return err;
#else
    printf(OPESCI_PAPI_MISSING);
#endif
}

int opesci_papi_start_counters(int numevents, int *events)
{
#if defined(OPESCI_HAVE_PAPI)
    if (PAPI_num_counters() < numevents) {
        printf("WARNING: More events specified that hardware counters available\n");
        printf(OPESCI_PAPI_WARN);
    }
    int err = PAPI_start_counters(events, numevents);
    if (err != PAPI_OK && omp_get_thread_num()==0) printf(OPESCI_PAPI_WARN);
#else
    printf(OPESCI_PAPI_MISSING);
#endif
}

int opesci_papi_name2event(char *name, int *event)
{
#if defined(OPESCI_HAVE_PAPI)
    return PAPI_event_name_to_code(name, event);
#else
    printf(OPESCI_PAPI_MISSING);
#endif
}

int opesci_papi_read_counters(int numevents, long long *counters)
{
#if defined(OPESCI_HAVE_PAPI)
    int err = PAPI_read_counters(counters, numevents);
    if (err != PAPI_OK && omp_get_thread_num()==0) printf(OPESCI_PAPI_WARN);
#else
    printf(OPESCI_PAPI_MISSING);
#endif
}

void opesci_flops(float *rtime, float *ptime, long long *flpins, float *mflops)
{
#if defined(OPESCI_HAVE_PAPI)
    int err = PAPI_flops(rtime, ptime, flpins, mflops);
    if (err != PAPI_OK && omp_get_thread_num()==0) printf(OPESCI_PAPI_WARN);
#else
    printf(OPESCI_PAPI_MISSING);
#endif
}
