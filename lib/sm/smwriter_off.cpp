/*
===============================================================================

  FILE:  smwriter_off.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2005  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    see corresponding header file

===============================================================================
*/
#include "smwriter_off.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "vec3fv.h"
#include "vec3iv.h"

void SMwriter_off::set_nverts(int nverts)
{
  this->nverts = nverts;
  vertex_buffer = (float*)realloc(vertex_buffer, sizeof(float)*3*nverts);
  vertex_buffer_alloc = nverts;
}

void SMwriter_off::set_nfaces(int nfaces)
{
  this->nfaces = nfaces;
  triangle_buffer = (int*)realloc(triangle_buffer, sizeof(int)*3*nfaces);
  triangle_buffer_alloc = nfaces;
}

bool SMwriter_off::open(FILE* file)
{
  if (file == 0)
  {
    return false;
  }

#ifdef _WIN32
  if (file == stdout)
  {
    if(_setmode( _fileno( stdout ), _O_TEXT ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdout to text (translated) mode\n");
    }
  }
#endif

  this->file = file;
  v_count = 0;
  f_count = 0;
  vertex_buffer_alloc = 1024;
  vertex_buffer = (float*)malloc(sizeof(float)*3*vertex_buffer_alloc);
  triangle_buffer_alloc = 2048;
  triangle_buffer = (int*)malloc(sizeof(int)*3*triangle_buffer_alloc);
  return (vertex_buffer != 0) && (triangle_buffer != 0);
}

void SMwriter_off::write_vertex(const float* v_pos_f)
{
  if (v_count == vertex_buffer_alloc)
  {
    vertex_buffer = (float*)realloc(vertex_buffer, sizeof(float)*3*vertex_buffer_alloc*2);
    vertex_buffer_alloc = vertex_buffer_alloc * 2;
  }
  VecCopy3fv(&(vertex_buffer[v_count*3]),v_pos_f);

  v_count++;
}

void SMwriter_off::write_triangle(const int* t_idx)
{
  if (f_count == triangle_buffer_alloc)
  {
    triangle_buffer = (int*)realloc(triangle_buffer, sizeof(int)*3*triangle_buffer_alloc*2);
    triangle_buffer_alloc = triangle_buffer_alloc * 2;
  }
  VecCopy3iv(&(triangle_buffer[f_count*3]),t_idx);

  f_count++;
}

void SMwriter_off::write_triangle(const int* t_idx, const bool* t_final)
{
  write_triangle(t_idx);
}

void SMwriter_off::write_finalized(int final_idx)
{
}

void SMwriter_off::close(bool close_file, bool update_header)
{
  // close of SMwriter_off
  int i;
  // write header
  fprintf(file, "OFF\012");
  fprintf(file, "%d %d 0\012",v_count,f_count);
  // write vertices
  for(i = 0; i < v_count; i++)
  {
    fprintf(file, "%.8g %.8g %.8g\012",vertex_buffer[3*i+0],vertex_buffer[3*i+1],vertex_buffer[3*i+2]);
  }
  // write triangles
  for(i = 0; i < f_count; i++)
  {
    fprintf(file, "3 %d %d %d\012",triangle_buffer[3*i+0],triangle_buffer[3*i+1],triangle_buffer[3*i+2]);
  }

  if (close_file && file /*&& file != stdout*/) fclose(file);
  file = 0;

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

void SMwriter_off::write_header()
{
  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      fprintf(file, "# %s\012",comments[i]);
    }
  }
  if (nverts != -1) fprintf(file, "# nverts %d\012",nverts);
  if (nfaces != -1) fprintf(file, "# nfaces %d\012",nfaces);
  if (bb_min_f) fprintf(file, "# bb_min %.8g %.8g %.8g\012",bb_min_f[0],bb_min_f[1],bb_min_f[2]);
  if (bb_max_f) fprintf(file, "# bb_max %.8g %.8g %.8g\012",bb_max_f[0],bb_max_f[1],bb_max_f[2]);
}

SMwriter_off::SMwriter_off()
{
  // init of SMwriter_off
  file = 0;
  vertex_buffer_alloc = 0;
  vertex_buffer = 0;
  triangle_buffer_alloc = 0;
  triangle_buffer = 0;
}

SMwriter_off::~SMwriter_off()
{
  // clean-up for SMwriter_off
  if (vertex_buffer) free(vertex_buffer);
  if (triangle_buffer) free(triangle_buffer);
}
