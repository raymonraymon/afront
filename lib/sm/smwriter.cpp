/*
===============================================================================

  FILE:  smwriter.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2003-07 martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    see corresponding header file

===============================================================================
*/
#include "smwriter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vec3fv.h"
#include "vec3iv.h"

void SMwriter::add_comment(const char* comment)
{
  if (comments == 0)
  {
    ncomments = 0;
    comments = (char**)malloc(sizeof(char*)*10);
    comments[9] = (char*)-1;
  }
  else if (comments[ncomments] == (char*)-1)
  {
    comments = (char**)realloc(comments,sizeof(char*)*ncomments*2);
    comments[ncomments*2-1] = (char*)-1;
  }
  comments[ncomments] = strdup(comment);
  ncomments++;
}

void SMwriter::set_nverts(int nverts)
{
  this->nverts = nverts;
}

void SMwriter::set_nfaces(int nfaces)
{
  this->nfaces = nfaces;
}

void SMwriter::set_boundingbox(const float* bb_min_f, const float* bb_max_f)
{
  if (this->bb_min_f == 0) this->bb_min_f = new float[3];
  if (this->bb_max_f == 0) this->bb_max_f = new float[3];
  VecCopy3fv(this->bb_min_f, bb_min_f);
  VecCopy3fv(this->bb_max_f, bb_max_f);
}

void SMwriter::close(bool close_file, bool update_header)
{
  fprintf(stderr, "WARNING: classes derived from SMwriter should implement their own close() function\n");
}

SMwriter::SMwriter()
{
  // init of SMwriter interface

  ncomments = 0;
  comments = 0;

  nverts = -1;
  nfaces = -1;

  v_count = -1;
  f_count = -1;

  bb_min_f = 0;
  bb_max_f = 0;

  need_pre_order = false;
}

SMwriter::~SMwriter()
{
  // clean-up for SMwriter interface
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
