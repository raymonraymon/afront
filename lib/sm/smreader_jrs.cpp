/*
===============================================================================

  FILE:  smreader_jrs.cpp
  
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
#include "smreader_jrs.h"

#include <stdlib.h>

#include "vec3fv.h"
#include "vec3iv.h"

#define MAX_LINE_SIZE 2048

bool SMreader_jrs::open(FILE* file_node, FILE* file_ele)
{
  char dummy[MAX_LINE_SIZE];
  int num_coords;
  int num_attrs;

  if (file_node == 0 || file_ele == 0)
  {
    return false;
  }

  skipped_lines = 0;

  // start processing the *.node file

  this->file_node = file_node;
  line_node = (char*)malloc(sizeof(char) * MAX_LINE_SIZE);
  if (fgets(line_node, sizeof(char) * MAX_LINE_SIZE, file_node) == 0)
  {
    free(line_node);
    line_node = 0;
    fprintf(stderr, "ERROR: running out of data when reading header of *.node file\n");
    return false;
  }

  // skip comments and empty lines of *.node file

  while (line_node[0] == '#' || line_node[1] == '#' || sscanf(line_node, "%s", dummy) < 1)
  {
    skipped_lines++;
    if (fgets(line_node, sizeof(char) * MAX_LINE_SIZE, file_node) == 0)
    {
      free(line_node);
      line_node = 0;
      fprintf(stderr, "ERROR: running out of data after skipping %d lines in *.node file\n", skipped_lines);
      return false;
    }
  }

  // read header of *.node file

  if (sscanf(line_node, "%d %d %d %s %s", &nverts, &num_coords, &num_attrs, dummy, dummy) != 4)
  {
    free(line_node);
    line_node = 0;
    fprintf(stderr, "ERROR: corrupt header in *.node file\n");
    return false;
  }
  
  ncoords = num_coords + num_attrs;

  // read the first line of actual data in *.node file

  if (fgets(line_node, sizeof(char) * MAX_LINE_SIZE, file_node) == 0)
  {
    free(line_node);
    line_node = 0;
    fprintf(stderr, "ERROR: running out of data in *.node file after reading 0 of %d points\n", nverts);
    return false;
  }

  // skip comments and empty lines of *.node file

  while (line_node[0] == '#' || line_node[1] == '#' || sscanf(line_node, "%s", dummy) < 1)
  {
    // we are skipping this line
    skipped_lines++;
    if (fgets(line_node, sizeof(char) * MAX_LINE_SIZE, file_node) == 0)
    {
      free(line_node);
      line_node = 0;
      fprintf(stderr, "ERROR: running out of data after skipping %d lines in *.node file\n", skipped_lines);
      return false;
    }
  }

  // start processing the ele file

  this->file_ele = file_ele;
  line_ele = (char*)malloc(sizeof(char) * MAX_LINE_SIZE);
  if (fgets(line_ele, sizeof(char) * MAX_LINE_SIZE, file_ele) == 0)
  {
    free(line_ele);
    line_ele = 0;
    fprintf(stderr, "ERROR: running out of data when reading header of *.ele file\n");
    return false;
  }

  // skip comments and empty lines of *.ele file

  while (line_ele[0] == '#' || line_ele[1] == '#' || sscanf(line_ele, "%s", dummy) < 1)
  {
    skipped_lines++;
    if (fgets(line_ele, sizeof(char) * MAX_LINE_SIZE, file_ele) == 0)
    {
      free(line_ele);
      line_ele = 0;
      fprintf(stderr, "ERROR: running out of data while skipping lines in *.ele file\n");
      return false;
    }
  }

  // read header of *.ele file

  if (sscanf(line_ele, "%d %s %s %s", &nfaces, dummy, dummy, dummy) != 3)
  {
    free(line_ele);
    line_ele = 0;
    fprintf(stderr, "ERROR: corrupt header in *.ele file\n");
    return false;
  }

  // read the first line of actual data in *.ele file

  if (fgets(line_ele, sizeof(char) * MAX_LINE_SIZE, file_ele) == 0)
  {
    free(line_ele);
    line_ele = 0;
    fprintf(stderr, "ERROR: running out of data in *.ele file after reading 0 of %d triangles\n", nfaces);
    return false;
  }

  // skip comments and empty lines after header in *.ele file

  while (line_ele[0] == '#' || line_ele[1] == '#' || sscanf(line_ele, "%s", dummy) < 1)
  {
    skipped_lines++;
    if (fgets(line_ele, sizeof(char) * MAX_LINE_SIZE, file_ele) == 0)
    {
      free(line_ele);
      line_ele = 0;
      fprintf(stderr, "ERROR: running out of data while skipping comments in *.ele file\n");
      return false;
    }
  }

  fprintf(stderr, "nverts: %d nfaces: %d (num_coords: %d num_attrs: %d)\n", nverts, nfaces, num_coords, num_attrs);

  v_count = 0;
  f_count = 0;

  v_idx = -1;
  final_idx = -1;
  v_pos_f[0] = v_pos_f[1] = v_pos_f[2] = 0.0f;
  t_idx[0] = t_idx[1] = t_idx[2] = -1;
  t_final[0] = t_final[1] = t_final[2] = false;

  return true;
}

bool SMreader_jrs::compute_bounding_box()
{
  int idx;
  // alloc bounding box
  if (bb_min_f == 0) bb_min_f = new float[3];
  if (bb_max_f == 0) bb_max_f = new float[3];
  // process first vertex
  if (ncoords >= 3)
  {
    sscanf(line_node, "%d %f %f %f", &idx, &(v_pos_f[0]), &(v_pos_f[1]), &(v_pos_f[2]));
  }
  else
  {
    sscanf(line_node, "%d %f %f", &idx, &(v_pos_f[0]), &(v_pos_f[1]));
    v_pos_f[2] = 0;
  }
  VecCopy3fv(bb_min_f, v_pos_f);
  VecCopy3fv(bb_max_f, v_pos_f);
  // process remaining vertices
  for (v_count = 1; v_count < nverts; v_count++)
  {
    if (fgets(line_node, sizeof(char) * MAX_LINE_SIZE, file_node) == 0)
    {
      free(line_node);
      line_node = 0;
      fprintf(stderr, "WARNING: running out of data in *.node file after reading %d of %d points\n", v_count, nverts);
      return true;
    }
    // read an element
    if (ncoords >= 3)
    {
      sscanf(line_node, "%d %f %f %f", &idx, &(v_pos_f[0]), &(v_pos_f[1]), &(v_pos_f[2]));
    }
    else
    {
      sscanf(line_node, "%d %f %f", &idx, &(v_pos_f[0]), &(v_pos_f[1]));
    }
    // update bounding box
    VecUpdateMinMax3fv(bb_min_f, bb_max_f, v_pos_f);
  }
  return true;
}

void SMreader_jrs::close(bool close_file)
{
  // close of SMreader_jrs
  if (skipped_lines) fprintf(stderr,"WARNING: skipped %d lines.\n",skipped_lines);
  skipped_lines = 0;
  if (close_file && file_node /*&& file_node != stdin*/) fclose(file_node);
  file_node = 0;
  if (close_file && file_ele /*&& file_ele != stdin*/) fclose(file_ele);
  file_ele = 0;
  if (line_node)
  {
    free(line_node);
    line_node = 0;
  }
  if (line_ele)
  {
    free(line_ele);
    line_ele = 0;
  }

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

SMevent SMreader_jrs::read_element()
{
  if (line_node == 0 && line_ele == 0)
  {
    return SM_EOF;
  }

  if (line_node && v_count < nverts)
  {
    v_idx = v_count;
    if (ncoords >= 3)
    {
      sscanf(line_node, "%d %f %f %f", &v_idx, &(v_pos_f[0]), &(v_pos_f[1]), &(v_pos_f[2]));
    }
    else
    {
      sscanf(line_node, "%d %f %f", &v_idx, &(v_pos_f[0]), &(v_pos_f[1]));
    }
    v_count++;
    if (fgets(line_node, sizeof(char) * MAX_LINE_SIZE, file_node) == 0)
    {
      if (v_count < nverts)
      {
        fprintf(stderr,"WARNING: early end of node file after %d of %d vertices\n", v_count, nverts);
      }
      free(line_node);
      line_node = 0;
    }
    return SM_VERTEX;
  }
  else if (line_ele && f_count < nfaces)
  {
    int idx;

    sscanf(line_ele, "%d %d %d %d", &idx, &(t_idx[0]), &(t_idx[1]), &(t_idx[2]));

    f_count++;

    if (fgets(line_ele, sizeof(char) * MAX_LINE_SIZE, file_ele) == 0)
    {
      if (f_count < nfaces)
      {
        fprintf(stderr,"WARNING: early end of ele file after %d of %d faces\n", f_count, nfaces);
      }
      free(line_ele);
      line_ele = 0;
    }
    return SM_TRIANGLE;
  }
  else
  {
    if (line_node)
    {
      free(line_node);
      line_node = 0;
    }
    if (line_ele)
    {
      free(line_ele);
      line_ele = 0;
    }
    return SM_EOF;
  }
}

SMevent SMreader_jrs::read_event()
{
  return read_element();
}

SMreader_jrs::SMreader_jrs()
{
  // init of SMreader_jrs
  ncoords = 0;
  skipped_lines = 0;
  file_node = 0;
  file_ele = 0;
  line_node = 0;
  line_ele = 0;
}

SMreader_jrs::~SMreader_jrs()
{
  // clean-up for SMreader_jrs
}
