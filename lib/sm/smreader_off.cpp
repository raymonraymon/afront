/*
===============================================================================

  FILE:  smreader_off.cpp
  
  CONTENTS:
  
    see corresponding header file
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005-07 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see corresponding header file
  
===============================================================================
*/
#include "smreader_off.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "vec3fv.h"
#include "vec3iv.h"

#define MAX_LINE_SIZE 1024

bool SMreader_off::open(FILE* file)
{
  if (file == 0)
  {
    return false;
  }

#ifdef _WIN32
  if (file == stdin)
  {
    if(_setmode( _fileno( stdin ), _O_TEXT ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdin to text (translated) mode\n");
    }
  }
#endif

  this->file = file;

  skipped_lines = 0;
  line = (char*)malloc(sizeof(char) * MAX_LINE_SIZE);
  if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
  {
    free(line);
    line = 0;
    return false;
  }

  post_order = false;

  bool has_off = false;

  if (strstr(line,"off") || strstr(line,"OFF"))
  {
    has_off = true;
    if (fgets(line, sizeof(char) * 256, file) == 0)
    {
      free(line);
      line = 0;
      return false;
    }
  }

  // look for header information
  if (sscanf(&(line[0]), "%d %d", &nverts, &nfaces) == 2)
  {
    fprintf(stderr, "nverts: %d nfaces: %d\n", nverts, nfaces);
  }
  else
  {
    fprintf(stderr, "ERROR: corrupt OFF file\n");
    return false;
  }

  // read next line
  if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
  {
    free(line);
    line = 0;
    return false;
  }
  else
  {
    v_count = 0;
  }

  fprintf(stderr, "nverts: %d nfaces: %d\n", nverts, nfaces);

  f_count = 0;

  v_idx = -1;
  final_idx = -1;
  v_pos_f[0] = v_pos_f[1] = v_pos_f[2] = 0.0f;
  t_idx[0] = t_idx[1] = t_idx[2] = -1;
  t_final[0] = t_final[1] = t_final[2] = false;

  return true;
}

bool SMreader_off::compute_bounding_box()
{
  // alloc bounding box
  if (bb_min_f == 0) bb_min_f = new float[3];
  if (bb_max_f == 0) bb_max_f = new float[3];
  // read first element
  if (sscanf(&(line[0]), "%f %f %f", &(v_pos_f[0]), &(v_pos_f[1]), &(v_pos_f[2])) != 3)
  {
    fprintf(stderr, "ERROR: corrupt OFF file\n");
    exit(1);
  }
  // read next line
  if (fgets(line, sizeof(char) * 256, file) == 0)
  {
    fprintf(stderr, "WARNING: early terminated OFF file\n");
    free(line);
    line = 0;
    return true;
  }
  // init bounding box
  VecCopy3fv(bb_min_f, v_pos_f);
  VecCopy3fv(bb_max_f, v_pos_f);
  // process remaining elements
  for (v_count = 1; v_count < nverts; v_count++)
  {
    // read next element
    if (sscanf(&(line[0]), "%f %f %f", &(v_pos_f[0]), &(v_pos_f[1]), &(v_pos_f[2])) != 3)
    {
      fprintf(stderr, "ERROR: corrupt OFF file\n");
      exit(1);
    }
    // read next line
    if (fgets(line, sizeof(char) * 256, file) == 0)
    {
      fprintf(stderr, "WARNING: early terminated OFF file\n");
      free(line);
      line = 0;
      return true;
    }
    // update bounding box
    VecUpdateMinMax3fv(bb_min_f, bb_max_f, v_pos_f);
  }
  return true;
}

void SMreader_off::close(bool close_file)
{
  // close of SMreader_off
  if (skipped_lines) fprintf(stderr,"WARNING: skipped %d lines.\n",skipped_lines);
  skipped_lines = 0;
  if (close_file && file /*&& file != stdin*/) fclose(file);
  file = 0;
  if (line)
  {
    free(line);
    line = 0;
  }
  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

SMevent SMreader_off::read_element()
{
  if (line == 0)
  {
    return SM_EOF;
  }

  if (v_count < nverts)
  {
    v_idx = v_count;
    while (sscanf(&(line[0]), "%f %f %f", &(v_pos_f[0]), &(v_pos_f[1]), &(v_pos_f[2])) != 3)
    {
      if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
      {
        fprintf(stderr,"WARNING: early end-of-file after %d of %d vertices\n", v_count, nverts);
        free(line);
        line = 0;
        return SM_EOF;
      }
      skipped_lines++;
    }
    v_count++;
    if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
    {
      if (v_count < nverts || f_count < nfaces)
      {
        fprintf(stderr,"WARNING: early end-of-file after %d of %d vertices\n", v_count, nverts);
      }
      free(line);
      line = 0;
    }
    return SM_VERTEX;
  }
  else if (f_count < nfaces)
  {
    int deg;
    while (true)
    {
      if (sscanf(&(line[0]), "%d %d %d %d", &(deg), &(t_idx[0]), &(t_idx[1]), &(t_idx[2])) == 4)
      {
        if (deg != 3) fprintf(stderr,"WARNING: polygon %d has degree of %d\n", f_count, deg);
        break;
      }
      else if (sscanf(&(line[0]), "%d %d %d", &(t_idx[0]), &(t_idx[1]), &(t_idx[2])) == 3)
      {
        break;
      }
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        if (f_count < nfaces)
        {
          fprintf(stderr,"WARNING: early end-of-file after %d of %d faces\n", f_count, nfaces);
        }
        free(line);
        line = 0;
        return SM_EOF;
      }
    }
    f_count++;
    if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
    {
      if (f_count < nfaces)
      {
        fprintf(stderr,"WARNING: early end-of-file after %d of %d faces\n", f_count, nfaces);
      }
      free(line);
      line = 0;
    }
    return SM_TRIANGLE;
  }
  free(line);
  line = 0;
  return SM_EOF;
}

SMevent SMreader_off::read_event()
{
  return read_element();
}

SMreader_off::SMreader_off()
{
  // init of SMreader_off
  file = 0;
  line = 0;
}

SMreader_off::~SMreader_off()
{
  // clean-up for SMreader_off
}
