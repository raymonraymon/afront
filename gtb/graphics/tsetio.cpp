
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


#include <gtb/gtb.hpp>
#ifndef WIN32
#include <gtb/error/error.hpp>
#include <gtb/memory/ptrs.h>
#include <gtb/replace/replace.h>
#endif

#include "tsetio.hpp"

#ifndef WIN32
#include <gtb/io/io.hpp>
#endif // WIN32

GTB_BEGIN_NAMESPACE

/*
 * A helper function that reads a face from a .obj file
 */
static void parse_obj_face(char* line, IndexedTriangleSet& triangles)
{
    const char* delim = " \t";
    char* token = strtok(line, delim);
    std::vector<unsigned> lflist; // local facelist
    int count=0;

    while (token != 0)
    {
        int v;

        if (sscanf(token, "%d", &v) < 1)
        {
            break;
        }

        if (v < 0) 
        {
#ifndef NO_EXCEPTIONS
            throw CErr("parse_obj_face: Cannot handle this type of .obj file");
#else
	    return;
#endif
        }

        --v; // vertex numbers in .obj file begin at 1
        lflist.push_back(v);

        token = strtok(0, delim);
        ++count;
    }

    if (count != 3)
#ifndef NO_EXCEPTIONS
      throw CErr("parse_obj_face: can handle triangles only");
#else
      return;
#endif

    triangles.insert_triangle(lflist[0], lflist[1], lflist[2], false, 0);
}       

void read_obj(const char* name, IndexedTriangleSet& triangles)
{
    afree<FILE*> f(xfopen(name, "r"), fclose);
    read_obj(f, triangles);
}

void read_obj(FILE* f, IndexedTriangleSet& triangles)
{
    while (!feof(f))
    {
        char line[1000];
        if (fgets(line, sizeof(line), f) == 0) break;


        if (strncmp(line, "vn ", 3) == 0)
        {
            double x,y,z;
            sscanf(line+2,"%lf %lf %lf", &x, &y, &z);
            triangles.insert_vertex_normal(Vector3(x,y,z));
        }
        else if (strncmp(line, "v ", 2) == 0)
        {
            double x,y,z;
            sscanf(line+1,"%lf %lf %lf", &x, &y, &z);
            triangles.insert_vertex(Point3(x,y,z));
        }
        else if (strncmp(line, "f ", 2) ==0)
        {
            parse_obj_face(line+2, triangles);
        }
    }
}

void write_obj(const char* name, const IndexedTriangleSet& triangles, int addidx)
{
    afree<FILE*> f(xfopen(name, "w"), fclose);

    {
        int n_vertices = triangles.num_vertices();
        for (int i = 0; i < n_vertices; ++i)
        {
            const Point3& v = triangles.vertex(i);
            fprintf(f, "v %f %f %f\n", v[0], v[1], v[2]);
        }
    }

    {
        int n_normals = triangles.num_vertex_normals();
        for (int i = 0; i < n_normals; ++i)
        {
            const Vector3& n = triangles.vertex_normal(i);
            fprintf(f, "vn %f %f %f\n", n[0], n[1], n[2]);
        }
    }

    {
        int n_faces = triangles.num_triangles();
        for (int i = 0; i < n_faces; ++i)
        {
            const IndexedTriangle& t = triangles.indexed_triangle(i);
            fprintf(f, "f %d %d %d\n", t.A()+addidx, t.B()+addidx, t.C()+addidx);
        }
    }
}

void read_triangleset(const char* name, IndexedTriangleSet& triangles)
{
    char ext[100];
    get_file_extension(name, ext, 100);

    if (stricmp(ext, "obj") == 0) read_obj(name, triangles);
    else if (stricmp(ext, "off") == 0)
    {
        afree<FILE*> f(xfopen(name, "rb"), fclose);
	    triangles.read_off(f);
    }
#ifndef NO_EXCEPTIONS
    else
      throw CErr("read_triangleset: Unknown file type");
#endif
}

GTB_END_NAMESPACE
