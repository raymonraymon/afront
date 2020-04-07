
/*

Copyright 2007 University of Utah


This file is part of Afront.

Afront is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Afront is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA

*/


// projectall.cpp : Defines the entry point for the console application.
//
// Demonstrates mlslib
//
// Takes as input a set of points and project them all on the surface
//

#include "stdafx.h"
using namespace mls;

typedef float areal;

char* me;
areal radius_factor;
typedef gtb::tsurfel_set<areal> surfel_set;
typedef gtb::tPoint3<areal> Point3;
typedef gtb::tVector3<areal> Vector3;

using gtb::CTimerOnce;

void usage()
{
    printf(
        "%s [-f] <src.obj> <dest.obj>\n"
        "  -f <factor> - radius factor, which defines smoothness default:%f\n",
        radius_factor,
        me);
    exit(-1);
}

void smooth_all(surfel_set& points, areal radius_factor, surfel_set& spoints);

int main(int argc, char** argv)
{
    --argc; ++argv;
    while (argc && (argv[0][0] == '-'))
    {
        switch(argv[0][1])
        {
        case 'f':
            radius_factor = atof(argv[1]);
            --argc; ++argv;
            break;
        default:
            usage();
        }
        --argc; ++argv;
    }

    if (argc != 2) usage();
    char* srcname = argv[0];
    char* destname = argv[1];

    surfel_set points;
    surfel_set spoints; // The smoothed points;
    gtb::read_points(srcname, points);
    smooth_all(points, radius_factor, spoints);
    write_points(destname, spoints);

	return 0;
}

void smooth_all(surfel_set& points, areal radius_factor, surfel_set& spoints)
{
    gtb::CTimerOnce __ttt("Projecting all samples", true);
    int n_failed = 0;

    GaussianWeightFunction<areal> theta;
    GaussianWeightFunction<areal> radius_wf;

    CProjection<areal> pobject(points, 8, &theta, &radius_wf, 2 /* degree */, 0);
    {
        CTimerOnce __t("Computing points radius");
        pobject.compute_points_radius();
    }

    pobject.set_radius_factor(radius_factor);

    ///
    // project all points
    //
    unsigned N = points.size();
    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& xi = points.vertex(i);
        Point3 r1;
        Vector3 n1;
        bool rc = pobject.PowellProject(xi, r1, n1);
        if (rc)
        {
            spoints.insert_vertex(r1, n1);
        }
        else ++n_failed;
        if (i%500 == 0)
        {
            printf("%d%%   (%d/%d)   \r", i*100/N, i, N);
        }
    }
    printf("                            \r");
    printf("Failed: %d\n", n_failed);
}

