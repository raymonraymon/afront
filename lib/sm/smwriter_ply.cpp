/*
===============================================================================

  FILE:  smwriter_ply.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2007  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    see corresponding header file

===============================================================================
*/
#include "smwriter_ply.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "vec3fv.h"
#include "vec3iv.h"

#include "ply.h"

/* vertex and face definitions for a polygonal object */

typedef struct Vertex {
  float x,y,z;
} Vertex;

typedef struct Face {
  unsigned char nverts;
  int *verts;
} Face;

static PlyProperty vertex_props[] = { /* list of property information for a vertex */
  {"x", Float32, Float32, offsetof(Vertex,x), 0, 0, 0, 0},
  {"y", Float32, Float32, offsetof(Vertex,y), 0, 0, 0, 0},
  {"z", Float32, Float32, offsetof(Vertex,z), 0, 0, 0, 0},
};

static PlyProperty face_props[] = { /* list of property information for a face */
  {"vertex_indices", Int32, Int32, offsetof(Face,verts), 1, Uint8, Uint8, offsetof(Face,nverts)},
};

void SMwriter_ply::set_nverts(int nverts)
{
  this->nverts = nverts;
  vertex_buffer = (float*)realloc(vertex_buffer, sizeof(float)*3*nverts);
  vertex_buffer_alloc = nverts;
}

void SMwriter_ply::set_nfaces(int nfaces)
{
  this->nfaces = nfaces;
  triangle_buffer = (int*)realloc(triangle_buffer, sizeof(int)*3*nfaces);
  triangle_buffer_alloc = nfaces;
}

bool SMwriter_ply::open(FILE* file)
{
  if (file == 0)
  {
    return false;
  }

#ifdef _WIN32
  if (file == stdout)
  {
    if(_setmode( _fileno( stdout ), _O_BINARY ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdout to binary (untranslated) mode\n");
    }
  }
#endif

  char* elem_names[2] = {"vertex","face"};

#ifdef _WIN32
  out_ply = write_ply(file, 2, elem_names, PLY_BINARY_LE);
#else
  out_ply = write_ply(file, 2, elem_names, PLY_BINARY_BE);
#endif

  if (out_ply == 0)
  {
    fprintf(stderr, "ERROR: cannot open output PLY file\n");
    return false;
  }

  v_count = 0;
  f_count = 0;
  vertex_buffer_alloc = 1024;
  vertex_buffer = (float*)malloc(sizeof(float)*3*vertex_buffer_alloc);
  triangle_buffer_alloc = 2048;
  triangle_buffer = (int*)malloc(sizeof(int)*3*triangle_buffer_alloc);
  return (vertex_buffer != 0) && (triangle_buffer != 0);
}

void SMwriter_ply::write_vertex(const float* v_pos_f)
{
  if (v_count == vertex_buffer_alloc)
  {
    vertex_buffer = (float*)realloc(vertex_buffer, sizeof(float)*3*vertex_buffer_alloc*2);
    vertex_buffer_alloc = vertex_buffer_alloc * 2;
  }
  VecCopy3fv(&(vertex_buffer[v_count*3]),v_pos_f);

  v_count++;
}

void SMwriter_ply::write_triangle(const int* t_idx)
{
  if (f_count == triangle_buffer_alloc)
  {
    triangle_buffer = (int*)realloc(triangle_buffer, sizeof(int)*3*triangle_buffer_alloc*2);
    triangle_buffer_alloc = triangle_buffer_alloc * 2;
  }
  VecCopy3iv(&(triangle_buffer[f_count*3]),t_idx);

  f_count++;
}

void SMwriter_ply::write_triangle(const int* t_idx, const bool* t_final)
{
  write_triangle(t_idx);
}

void SMwriter_ply::write_finalized(int final_idx)
{
}

void SMwriter_ply::close(bool close_file, bool update_header)
{
  // close of SMwriter_ply
  int i;
  // write header
  element_layout_ply(((PlyFile*)out_ply), "vertex", v_count, 3, vertex_props);
  element_layout_ply(((PlyFile*)out_ply), "face", f_count, 1, face_props);
  header_complete_ply(((PlyFile*)out_ply));
  // write vertices
  put_element_setup_ply(((PlyFile*)out_ply), "vertex");
  for(i = 0; i < v_count; i++)
  {
    put_element_ply(((PlyFile*)out_ply),(void*)&(vertex_buffer[3*i]));
  }
  // write triangles
  Face face;
  face.nverts = 3;
  put_element_setup_ply(((PlyFile*)out_ply), "face");
  for(i = 0; i < f_count; i++)
  {
    face.verts = &(triangle_buffer[3*i]);
    put_element_ply(((PlyFile*)out_ply),(void*)&face);
  }

  if (close_file && out_ply && ((PlyFile*)out_ply)->fp /*&& ((PlyFile*)out_ply)->fp != stdout*/) fclose(((PlyFile*)out_ply)->fp);
  if (out_ply) free_ply(((PlyFile*)out_ply));
  out_ply = 0;

  vertex_buffer_alloc = 0;
  if (vertex_buffer) free(vertex_buffer);
  vertex_buffer = 0;

  triangle_buffer_alloc = 0;
  if (triangle_buffer) free(triangle_buffer);
  triangle_buffer = 0;

  // close of SMwriter interface
  if (nverts != -1 && nverts != v_count)  fprintf(stderr,"WARNING: set nverts %d but v_count %d\n",nverts,v_count);
  if (nfaces != -1 && nfaces != f_count)  fprintf(stderr,"WARNING: set nfaces %d but f_count %d\n",nfaces,f_count);
  nverts = v_count;
  nfaces = f_count;
  v_count = -1;
  f_count = -1;
}

void SMwriter_ply::write_header()
{
  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      fprintf(((PlyFile*)out_ply)->fp, "# %s\012",comments[i]);
    }
  }
  if (nverts != -1) fprintf(((PlyFile*)out_ply)->fp, "# nverts %d\012",nverts);
  if (nfaces != -1) fprintf(((PlyFile*)out_ply)->fp, "# nfaces %d\012",nfaces);
  if (bb_min_f) fprintf(((PlyFile*)out_ply)->fp, "# bb_min %.8g %.8g %.8g\012",bb_min_f[0],bb_min_f[1],bb_min_f[2]);
  if (bb_max_f) fprintf(((PlyFile*)out_ply)->fp, "# bb_max %.8g %.8g %.8g\012",bb_max_f[0],bb_max_f[1],bb_max_f[2]);
}

SMwriter_ply::SMwriter_ply()
{
  // init of smwriter_ply interface
  out_ply = 0;
  vertex_buffer_alloc = 0;
  vertex_buffer = 0;
  triangle_buffer_alloc = 0;
  triangle_buffer = 0;
}

SMwriter_ply::~SMwriter_ply()
{
  // clean-up  of smwriter_ply interface
  if (out_ply) free_ply (((PlyFile*)out_ply));
  if (vertex_buffer) free(vertex_buffer);
  if (triangle_buffer) free(triangle_buffer);
}
