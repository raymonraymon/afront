/*
===============================================================================

  FILE:  smreader_ply.cpp
  
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
#include "smreader_ply.h"

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

bool SMreader_ply::open(FILE* file)
{
  int i;
  int elem_count;
  char *elem_name;
  
  if (file == 0)
  {
    return false;
  }

#ifdef _WIN32
  if (file == stdin)
  {
    if(_setmode( _fileno( stdin ), _O_BINARY ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdin to binary (untranslated) mode\n");
    }
  }
#endif

  in_ply = read_ply(file);

  if (in_ply == 0)
  {
    fprintf(stderr, "FATAL ERROR: input PLY file is corrupt\n");
    return false;
  }

  nverts = 0;
  nfaces = 0;

  v_count = 0;
  f_count = 0;

  for (i = 0; i < ((PlyFile*)in_ply)->num_elem_types; i++)
  {
    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply(((PlyFile*)in_ply), i, &elem_count);

    if (equal_strings ("vertex", elem_name))
    {
      nverts = elem_count;
      velem = i;
    }
    else if (equal_strings ("face", elem_name))
    {
      nfaces = elem_count;
      felem = i;
    }
    else
    {
      fprintf(stderr, "WARNING: unknown element '%s' in PLY file\n", elem_name);
    }
  }
  
  if (nverts == 0)
  {
    fprintf(stderr, "WARNING: no vertices in PLY file\n");
  }
  else
  {
    /* prepare to read the vertex elements */
    setup_element_read_ply(((PlyFile*)in_ply), velem, &elem_count);
    /* set up for getting vertex elements */
    setup_property_ply (((PlyFile*)in_ply), &vertex_props[0]);
    setup_property_ply (((PlyFile*)in_ply), &vertex_props[1]);
    setup_property_ply (((PlyFile*)in_ply), &vertex_props[2]);
  }

  if (nfaces == 0)
  {
    fprintf(stderr, "WARNING: no triangles in PLY file\n");
  }

  v_idx = -1;
  final_idx = -1;
  v_pos_f[0] = v_pos_f[1] = v_pos_f[2] = 0.0f;
  t_idx[0] = t_idx[1] = t_idx[2] = -1;
  t_final[0] = t_final[1] = t_final[2] = false;

  return true;
}

bool SMreader_ply::compute_bounding_box()
{
  // alloc bounding box
  if (bb_min_f == 0) bb_min_f = new float[3];
  if (bb_max_f == 0) bb_max_f = new float[3];
  // read first element
  get_element_ply (((PlyFile*)in_ply), (void *)v_pos_f);
  // init bounding box
  VecCopy3fv(bb_min_f, v_pos_f);
  VecCopy3fv(bb_max_f, v_pos_f);
  // process remaining elements
  for (v_count = 1; v_count < nverts; v_count++)
  {
    // read an element
    get_element_ply (((PlyFile*)in_ply), (void *)v_pos_f);
    // update bounding box
    VecUpdateMinMax3fv(bb_min_f, bb_max_f, v_pos_f);
  }
  return true;
}

void SMreader_ply::close(bool close_file)
{
  // close of SMreader_ply
	if (close_file && in_ply && ((PlyFile*)in_ply)->fp /*&& ((PlyFile*)in_ply)->fp != stdin*/) fclose(((PlyFile*)in_ply)->fp);
  if (in_ply) free_ply(((PlyFile*)in_ply));

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

SMevent SMreader_ply::read_element()
{
  int elem_count;

  if (v_count < nverts)
  {
    v_idx = v_count;
    get_element_ply(((PlyFile*)in_ply), (void *)v_pos_f);
    v_count++;
    return SM_VERTEX;
  }
  else if (f_count < nfaces)
  {
    if (f_count == 0)
    {
      /* prepare to read the face elements */
      setup_element_read_ply(((PlyFile*)in_ply), felem, &elem_count);
      /* set up for getting face elements */
      setup_property_ply (((PlyFile*)in_ply), &face_props[0]);
    }
    Face* face = (Face *) malloc (sizeof (Face));
    get_element_ply (((PlyFile*)in_ply), (void *) face);
    VecCopy3iv(t_idx, face->verts);
    free(face->verts);
    free(face);
    f_count++;
    return SM_TRIANGLE;
  }
  else
  {
    return SM_EOF;
  }
}

SMevent SMreader_ply::read_event()
{
  return read_element();
}

SMreader_ply::SMreader_ply()
{
  // init of SMreader_ply
  in_ply = 0;
  velem = -1;
  felem = -1;
}

SMreader_ply::~SMreader_ply()
{
  // clean-up for SMreader_ply
}
