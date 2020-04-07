/*
===============================================================================

  FILE:  smreader.cpp
  
  CONTENTS:
  
    see corresponding header file
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003-07  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see corresponding header file
  
===============================================================================
*/
#include "smreader.h"

#include <stdio.h>
#include <stdlib.h>

#include "vec3fv.h"

bool SMreader::compute_bounding_box()
{
  SMevent event;

  if (bb_min_f == 0) bb_min_f = new float[3];
  if (bb_max_f == 0) bb_max_f = new float[3];

  while ((event = read_element()) > SM_EOF)
  {
    if (event == SM_VERTEX)
    {
      VecCopy3fv(bb_min_f, v_pos_f);
      VecCopy3fv(bb_max_f, v_pos_f);
      break;
    }
  }
  while ((event = read_element()) > SM_EOF)
  {
    if (event == SM_VERTEX)
    {
      VecUpdateMinMax3fv(bb_min_f, bb_max_f, v_pos_f);
    }
  }

  return true;
}

void SMreader::close(bool close_file)
{
  fprintf(stderr, "WARNING: classes derived from SMreader should implement their own close() function\n");
}

SMreader::SMreader()
{
  // init of SMreader interface

  v_idx = -1;
  v_pos_f[0] = v_pos_f[1] = v_pos_f[2] = 0;

  t_idx[0] = t_idx[1] = t_idx[2] = -1;
  t_final[0] = t_final[1] = t_final[2] = false;

  final_idx = -1;

  ncomments = 0;
  comments = 0;

  nfaces = -1;
  nverts = -1;

  f_count = -1;
  v_count = -1;

  bb_min_f = 0;
  bb_max_f = 0;

  post_order = false;
}

SMreader::~SMreader()
{
  // clean-up for SMreader interface
  if (v_count != -1)
  {
    close(); // user must have forgotten to close the mesh
  }
  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      free(comments[i]);
    }
    free(comments);
  }
  if (bb_min_f) delete [] bb_min_f;
  if (bb_max_f) delete [] bb_max_f;
}
