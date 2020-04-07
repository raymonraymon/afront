/*
===============================================================================

  FILE:  smwriter_sma.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2003  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    see corresponding header file

===============================================================================
*/
#include "smwriter_sma.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "vec3fv.h"
#include "vec3iv.h"

bool SMwriter_sma::open(FILE* file)
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
  return true;
}

void SMwriter_sma::write_vertex(const float* v_pos_f)
{
  if (v_count + f_count == 0) write_header();

  fprintf(file, "v %.8g %.8g %.8g\012",v_pos_f[0],v_pos_f[1],v_pos_f[2]);
  v_count++;
}

void SMwriter_sma::write_triangle(const int* t_idx)
{
  if (v_count + f_count == 0) write_header();

  fprintf(file, "f %d %d %d\012",t_idx[0]+1,t_idx[1]+1,t_idx[2]+1);
  f_count++;
}

void SMwriter_sma::write_triangle(const int* t_idx, const bool* t_final)
{
  if (v_count + f_count == 0) write_header();

  fprintf(file, "f %d %d %d\012",(t_final[0] ? t_idx[0]-v_count : t_idx[0]+1),(t_final[1] ? t_idx[1]-v_count : t_idx[1]+1), (t_final[2] ? t_idx[2]-v_count : t_idx[2]+1));
  f_count++;
}

void SMwriter_sma::write_finalized(int final_idx)
{
  if (final_idx < 0)
  {
    fprintf(file, "x %d\012",final_idx);
  }
  else
  {
    fprintf(file, "x %d\012",final_idx+1);
  }
}

void SMwriter_sma::close(bool close_file, bool update)
{
  // close of SMwriter_sma
  if (nverts != v_count || nfaces != f_count)
  {
	  if (update /*&& file != stdout*/) update_header();
  }

  if (close_file && file /*&& file != stdout*/) fclose(file);
  file = 0;

  // close of SMwriter interface
  if (nverts != -1 && nverts != v_count)  fprintf(stderr,"WARNING: set nverts %d but v_count %d\n",nverts,v_count);
  if (nfaces != -1 && nfaces != f_count)  fprintf(stderr,"WARNING: set nfaces %d but f_count %d\n",nfaces,f_count);
  nverts = v_count;
  nfaces = f_count;
  v_count = -1;
  f_count = -1;
}

void SMwriter_sma::write_header()
{
	if (nverts != -1/* || file != stdout*/) fprintf(file, "# nverts %-12d\012",nverts);
	if (nfaces != -1/* || file != stdout*/) fprintf(file, "# nfaces %-12d\012",nfaces);
  if (bb_min_f) fprintf(file, "# bb_min %.8g %.8g %.8g\012",bb_min_f[0],bb_min_f[1],bb_min_f[2]);
  if (bb_max_f) fprintf(file, "# bb_max %.8g %.8g %.8g\012",bb_max_f[0],bb_max_f[1],bb_max_f[2]);
  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      fprintf(file, "# %s\012",comments[i]);
    }
  }
}

void SMwriter_sma::update_header()
{
	/*
  if (file == stdout)
  {
    fprintf(stderr,"WARNING: cannot update header when writing to stdout\n");
    return;
  }
  if (fseek(file, 0L, SEEK_SET) != 0)
  {
    fprintf(stderr,"WARNING: cannot update header because fseek failed\n");
    return;
  }
  fprintf(file, "# nverts %-12d\012",v_count);
  fprintf(file, "# nfaces %-12d\012",f_count);
  if (fseek(file, 0L, SEEK_END) != 0)
  {
    fprintf(stderr,"WARNING: header was updated but final fseek failed\n");
  }
	*/
}

SMwriter_sma::SMwriter_sma()
{
  // init of SMwriter_sma
  file = 0;
}

SMwriter_sma::~SMwriter_sma()
{
  // clean-up for SMwriter_sma
}
