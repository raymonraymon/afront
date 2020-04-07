/*
===============================================================================

  FILE:  smreadwithsecondderivatives.cpp
  
  CONTENTS:
  
    see corresponding header file

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2007 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see corresponding header file
  
===============================================================================
*/
#include "smreadwithsecondderivatives.h"

#include <stdio.h>
#include "hash_map.h"
#include "deque.h"
#include <assert.h>

#include "vec3fv.h"

#define PRINT_CONTROL_OUTPUT
#undef PRINT_CONTROL_OUTPUT

struct SMtriangle;

typedef struct SMvertex
{
  SMvertex* buffer_next;    // used for efficient memory management
  float v[3];
  float n[3];
  float dz[3];
  float d2z[3];
  float nDx[3];
  float nDy[3];
  int index;
  short list_size;
  short list_alloc;
  SMtriangle** list;
} SMvertex;

typedef struct SMtriangle
{
  SMtriangle* buffer_next;  // used for efficient memory management
  SMvertex* vertices[3];
  char finalized_vertices;
} SMtriangle;

typedef hash_map<int, SMvertex*> my_vertex_hash;
typedef deque<SMtriangle*> my_triangle_queue;

static my_vertex_hash* vertex_hash;
static my_triangle_queue* triangle_queue;

// efficient memory allocation

static SMvertex* vertex_buffer_next;
static int vertex_buffer_increment;
static int vertex_buffer_number_chunks;
static SMvertex* vertex_buffer_chunks[16];

#ifdef PRINT_CONTROL_OUTPUT
static int vertex_buffer_size;
static int vertex_buffer_maxsize;
#endif // PRINT_CONTROL_OUTPUT

static void initVertexBuffer(int increment)
{
  vertex_buffer_next = 0;
  vertex_buffer_increment = increment;
  vertex_buffer_number_chunks = 0;
#ifdef PRINT_CONTROL_OUTPUT
  vertex_buffer_size = 0;
  vertex_buffer_maxsize = 0;
#endif // PRINT_CONTROL_OUTPUT
}

static SMvertex* allocVertex()
{
  if (vertex_buffer_next == 0)
  {
    vertex_buffer_next = (SMvertex*)malloc(sizeof(SMvertex)*vertex_buffer_increment);
    if (vertex_buffer_next == 0)
    {
      fprintf(stderr,"FATAL ERROR: malloc for vertex buffer failed\n");
      return 0;
    }
    else
    {
      vertex_buffer_chunks[vertex_buffer_number_chunks] = vertex_buffer_next;
      vertex_buffer_number_chunks++;
    }
    for (int i = 0; i < vertex_buffer_increment; i++)
    {
      vertex_buffer_next[i].buffer_next = &(vertex_buffer_next[i+1]);
      vertex_buffer_next[i].list = 0;
    }
    vertex_buffer_next[vertex_buffer_increment-1].buffer_next = 0;
    vertex_buffer_increment *= 2;
  }
  // get pointer to next available vertex
  SMvertex* vertex = vertex_buffer_next;
  vertex_buffer_next = vertex->buffer_next;

  // initialize data fields
  if (vertex->list == 0)
  {
    vertex->list = (SMtriangle**)malloc(sizeof(SMtriangle*)*10);
    vertex->list_alloc = 10;
  }
  vertex->list_size = 0;

#ifdef PRINT_CONTROL_OUTPUT
  vertex_buffer_size++;
  if (vertex_buffer_size > vertex_buffer_maxsize) vertex_buffer_maxsize = vertex_buffer_size;
#endif // PRINT_CONTROL_OUTPUT

  return vertex;
}

static void deallocVertex(SMvertex* vertex)
{
  vertex->buffer_next = vertex_buffer_next;
  vertex_buffer_next = vertex;
#ifdef PRINT_CONTROL_OUTPUT
  vertex_buffer_size--;
#endif // PRINT_CONTROL_OUTPUT
}

static void destroyVertexBuffer()
{
  for (int i = 0; i < vertex_buffer_number_chunks; i++)
  {
    free(vertex_buffer_chunks[i]);
  }
  vertex_buffer_number_chunks = 0;
  vertex_buffer_next = 0;
}

static SMtriangle* triangle_buffer_next;
static int triangle_buffer_increment;
static int triangle_buffer_number_chunks;
static SMtriangle* triangle_buffer_chunks[16];

#ifdef PRINT_CONTROL_OUTPUT
static int triangle_buffer_size;
static int triangle_buffer_maxsize;
#endif // PRINT_CONTROL_OUTPUT

static void initTriangleBuffer(int increment)
{
  triangle_buffer_next = 0;
  triangle_buffer_increment = increment;
  triangle_buffer_number_chunks = 0;
#ifdef PRINT_CONTROL_OUTPUT
  triangle_buffer_size = 0;
  triangle_buffer_maxsize = 0;
#endif // PRINT_CONTROL_OUTPUT
}

static SMtriangle* allocTriangle()
{
  if (triangle_buffer_next == 0)
  {
    triangle_buffer_next = (SMtriangle*)malloc(sizeof(SMtriangle)*triangle_buffer_increment);
    if (triangle_buffer_next == 0)
    {
      fprintf(stderr,"FATAL ERROR: malloc for triangle buffer failed\n");
      return 0;
    }
    else
    {
      triangle_buffer_chunks[triangle_buffer_number_chunks] = triangle_buffer_next;
      triangle_buffer_number_chunks++;
    }
    for (int i = 0; i < triangle_buffer_increment; i++)
    {
      triangle_buffer_next[i].buffer_next = &(triangle_buffer_next[i+1]);
    }
    triangle_buffer_next[triangle_buffer_increment-1].buffer_next = 0;
    triangle_buffer_increment *= 2;
  }
  // get pointer to next available triangle
  SMtriangle* triangle = triangle_buffer_next;
  triangle_buffer_next = triangle->buffer_next;
 
  // initialize data fields
  triangle->finalized_vertices = 0;

#ifdef PRINT_CONTROL_OUTPUT
  triangle_buffer_size++;
  if (triangle_buffer_size > triangle_buffer_maxsize) triangle_buffer_maxsize = triangle_buffer_size;
#endif // PRINT_CONTROL_OUTPUT

  return triangle;
}

static void deallocTriangle(SMtriangle* triangle)
{
  triangle->buffer_next = triangle_buffer_next;
  triangle_buffer_next = triangle;
#ifdef PRINT_CONTROL_OUTPUT
  triangle_buffer_size--;
#endif // PRINT_CONTROL_OUTPUT
}

static void destroyTriangleBuffer()
{
  for (int i = 0; i < triangle_buffer_number_chunks; i++)
  {
    free(triangle_buffer_chunks[i]);
  }
  triangle_buffer_number_chunks = 0;
  triangle_buffer_next = 0;
}

static void contribute_triangle_derivatives(SMvertex** vertices)
{
  float u[3];
  float v[3];
  float c[3];
  // add normal to vertices
  for (int i = 0; i < 3; i++)
  {
    SMvertex *vertex = vertices[i];
    SMvertex *v2 = vertices[(i+2) % 3];
    SMvertex *v3 = vertices[(i+1) % 3];

    // Calculate.
    //Vector3 ux(v2->v[0] - vertex->v[0],
    //       v2->v[1] - vertex->v[1],
    //       v2->dzdx - vertex->dzdx);
    //Vector3 vx(v3->v[0] - vertex->v[0],
    //       v3->v[1] - vertex->v[1],
    //       v3->dzdx - vertex->dzdx);
    //nDx += Cross(ux, vx);
    u[0] = v2->v[0] - vertex->v[0];
    u[1] = v2->v[1] - vertex->v[1];
    u[2] = v2->dz[0] - vertex->dz[0];
    
    v[0] = v3->v[0] - vertex->v[0];
    v[1] = v3->v[1] - vertex->v[1];
    v[2] = v3->dz[0] - vertex->dz[0];
      
    VecCrossProd3fv(c, u, v);
    VecSelfAdd3fv(vertex->nDx, c);
      
    //Vector3 uy(v2->v[0] - vertex->v[0],
    //       v2->v[1] - vertex->v[1],
    //       v2->dzdy - vertex->dzdy);
    //Vector3 vy(v3->v[0] - vertex->v[0],
    //       v3->v[1] - vertex->v[1],
    //       v3->dzdy - vertex->dzdy);
    //nDy += Cross(uy, vy);
    u[0] = v2->v[0] - vertex->v[0];
    u[1] = v2->v[1] - vertex->v[1];
    u[2] = v2->dz[1] - vertex->dz[1];
    
    v[0] = v3->v[0] - vertex->v[0];
    v[1] = v3->v[1] - vertex->v[1];
    v[2] = v3->dz[1] - vertex->dz[1];
    
    VecCrossProd3fv(c, u, v);
    VecSelfAdd3fv(vertex->nDy, c);
  }
}

static void finalize_vertex(SMvertex* vertex)
{
  // compute second derivatives
  vertex->d2z[0] = -(vertex->nDx[0] / vertex->nDx[2]); // d2z/dx2
  vertex->d2z[1] = -(vertex->nDy[1] / vertex->nDy[2]); // d2z/dy2
  vertex->d2z[2] = 0.5f * (-vertex->nDx[1] / vertex->nDx[2] - vertex->nDy[0] / vertex->nDy[2]); // d2z/dxdy
  // loop over the triangles of this vertex
  SMtriangle* triangle;
  for (int i = 0; i < vertex->list_size; i++)
  {
    // get the next triangle
    triangle = vertex->list[i];
    // increment that triangles finalized vertex counter
    triangle->finalized_vertices++;
    // are all vertices of this triangle finalized
    if (triangle->finalized_vertices == 3)
    {
      triangle_queue->push_back(triangle);
    }
  }
}

bool SMreadWithSecondDerivatives::open(SMreader* smreader)
{
  if (smreader == 0)
  {
    return false;
  }
  SMreadWithVertexNormals* smreadwithvertexnormals = new SMreadWithVertexNormals();
  if (smreadwithvertexnormals->open(smreader) == false)
  {
    return false;
  }
  return open(smreadwithvertexnormals);
}

bool SMreadWithSecondDerivatives::open(SMreadWithVertexNormals* smreadwithvertexnormals)
{
  if (smreadwithvertexnormals == 0)
  {
    return false;
  }
  this->smreadwithvertexnormals = smreadwithvertexnormals;

  nverts = smreadwithvertexnormals->nverts;
  nfaces = smreadwithvertexnormals->nfaces;

  v_count = 0;
  f_count = 0;

  bb_min_f = smreadwithvertexnormals->bb_min_f;
  bb_max_f = smreadwithvertexnormals->bb_max_f;

  vertex_hash = new my_vertex_hash;
  triangle_queue = new my_triangle_queue;

  initVertexBuffer(1024);
  initTriangleBuffer(2048);

  return true;
}

void SMreadWithSecondDerivatives::close(bool close_file)
{
  // close of SMreadWithSecondDerivatives
  smreadwithvertexnormals->close(close_file);

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"vertex_buffer_maxsize (size) %d (%d) triangle_buffer_maxsize (size) t %d (%d)\n",vertex_buffer_maxsize,vertex_buffer_size,triangle_buffer_maxsize,triangle_buffer_size);
#endif

  delete vertex_hash;
  delete triangle_queue;

  destroyVertexBuffer();
  destroyTriangleBuffer();

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

int SMreadWithSecondDerivatives::read_input()
{
  int i;
  SMvertex* vertex;
  SMtriangle* triangle;
  my_vertex_hash::iterator hash_elements[3];

  SMevent event = smreadwithvertexnormals->read_element();
    
  if (event == SM_TRIANGLE)
  {
    // create triangle
    triangle = allocTriangle();
    // lookup the triangle's vertices
    for (i = 0; i < 3; i++)
    {
      // look for vertex in the hash
      hash_elements[i] = vertex_hash->find(smreadwithvertexnormals->t_idx[i]);
      // did we find the vertex in the hash  
      if (hash_elements[i] == vertex_hash->end())
      {
        // fatal error
        fprintf(stderr, "FATAL ERROR: vertex %d not in hash. corrupt smreadwithvertexnormals output.\n", smreadwithvertexnormals->t_idx[i]);
        return -1;
      }
      else
      {
        // use vertex found in hash
        vertex = (*hash_elements[i]).second;
      }
      // add vertex to triangle
      triangle->vertices[i] = vertex;
      // add triangle to vertex
      if (vertex->list_size == vertex->list_alloc)
      {
        vertex->list_alloc += 4;
        vertex->list = (SMtriangle**)realloc(vertex->list,vertex->list_alloc*sizeof(SMtriangle*));
      }
      vertex->list[vertex->list_size++] = triangle;
    }
    // compute and add its normal to its vertices
    contribute_triangle_derivatives(triangle->vertices);
    // check for finalization of the triangle's vertices
    for (i = 0; i < 3; i++)
    {
      // look for vertex in the hash
      if (smreadwithvertexnormals->t_final[i])
      {
        // remove vertex from hash
        vertex_hash->erase(hash_elements[i]);
        // finalize vertex
        finalize_vertex(triangle->vertices[i]);
      }
    }
  }
  else if (event == SM_VERTEX)
  {
    // create the vertex ...
    vertex = allocVertex();
    // ... insert it into hash
    vertex_hash->insert(my_vertex_hash::value_type(smreadwithvertexnormals->v_idx, vertex));
    // ... copy coordinates ...
    VecCopy3fv(vertex->v, smreadwithvertexnormals->v_pos_f);
    // ... copy normalized normal ...
    VecNormalize3fv(vertex->n, smreadwithvertexnormals->v_nor_f);
    // ... compute first derivative ...
    vertex->dz[0] = -vertex->n[0] / vertex->n[2];
    vertex->dz[1] = -vertex->n[1] / vertex->n[2];
    // ... zero the accumulating Dx and Dy
    VecZero3fv(vertex->nDx);
    VecZero3fv(vertex->nDy);
    // ... init triangle list
    vertex->list_size = 0;
  }
  else if (event == SM_FINALIZED)
  {
    // look for vertex in the vertex hash
    hash_elements[0] = vertex_hash->find(smreadwithvertexnormals->v_idx);
    // did we find the vertex in the hash  
    if (hash_elements[0] == vertex_hash->end())
    {
      // fatal error
      fprintf(stderr, "FATAL ERROR: finalized vertex %d not in hash. corrupt streaming mesh.\n", smreadwithvertexnormals->v_idx);
      return -1;
    }
    else
    {
      // use vertex found in hash
      vertex = vertex = (*hash_elements[0]).second;
    }
    // we finalize it ...
    finalize_vertex(vertex);
    // ... and remove it from hash
    vertex_hash->erase(hash_elements[0]);
  }
  else if (event == SM_EOF)
  {
    if (vertex_hash->size())
    {
      // get some hash element ... 
      hash_elements[0] = vertex_hash->begin();
      // ... or rather its unfinalized vertex ...
      vertex = vertex = (*hash_elements[0]).second;
      // ... and finalize it ...
      finalize_vertex(vertex);
      // ... and remove it from hash
      vertex_hash->erase(hash_elements[0]);
    }
    else
    {
      return 0;
    }
  }
  return 1;
}

int SMreadWithSecondDerivatives::read_triangle()
{
  int i;
  SMvertex* vertex;
  SMtriangle* triangle;

  // read from input until we have a triangle for output (or fail)

  while (!triangle_queue->size())
  {
    int ok = read_input();
    if (ok <= 0)
    {
      return ok;
    }
  }

  // get next output triangle

  triangle = triangle_queue->front();
  triangle_queue->pop_front();

  // make sure it has everything

  assert(triangle->finalized_vertices == 3);
  
  // now we have a triangle

  have_triangle = 1;

  // look for new or finalized vertices and fill in indices, finalization, and positions

  for (i = 0; i < 3; i++)
  {
    vertex = triangle->vertices[i];
    // is this vertex a new vertex?
    if (vertex->index == -1)
    {
      vertex->index = v_count + have_new;
      new_vertices[have_new++] = i;
    }
    assert(vertex->index >= 0);
    // is this vertex to be finalized?
    if (vertex->list_size == 1)
    {
      t_final[i] = true;
      finalized_vertices[have_finalized++] = i;
      deallocVertex(vertex); // the memory can still be used until the next vertex is allocated.
    }
    else
    {
      vertex->list_size--;
      t_final[i] = false;
    }
    // fill index array
    t_idx[i] = vertex->index;
    // set additional position array
    t_pos_f[i] = vertex->v;
    // set additional normal array
    t_nor_f[i] = vertex->n;
  }

  deallocTriangle(triangle);

  return 1;
}

SMevent SMreadWithSecondDerivatives::read_element()
{
  if ((have_new + have_triangle) == 0)
  {
    next_new = 0;
    have_finalized = next_finalized = 0;
    
    int ok = read_triangle();
    if (ok <= 0)
    {
      return (SMevent)ok;
    }
  }
  if (have_new)
  {
    v_idx = t_idx[new_vertices[next_new]];
    VecCopy3fv(v_pos_f, t_pos_f[new_vertices[next_new]]);
    VecCopy3fv(v_nor_f, t_nor_f[new_vertices[next_new]]);
    VecCopy3fv(v_dz_f, t_dz_f[new_vertices[next_new]]);
    VecCopy3fv(v_d2z_f, t_d2z_f[new_vertices[next_new]]);
    have_new--; next_new++;
    v_count++;
    return SM_VERTEX;
  }
  else
  {
    have_triangle = 0;
    f_count++;
    return SM_TRIANGLE;
  }
}

SMevent SMreadWithSecondDerivatives::read_event()
{
  if ((have_new + have_triangle + have_finalized) == 0)
  {
    return read_element();
  }
  if (have_new)
  {
    v_idx = t_idx[new_vertices[next_new]];
    VecCopy3fv(v_pos_f, t_pos_f[new_vertices[next_new]]);
    VecCopy3fv(v_nor_f, t_nor_f[new_vertices[next_new]]);
    VecCopy3fv(v_dz_f, t_dz_f[new_vertices[next_new]]);
    VecCopy3fv(v_d2z_f, t_d2z_f[new_vertices[next_new]]);
    have_new--; next_new++;
    v_count++;
    return SM_VERTEX;
  }
  else if (have_triangle)
  {
    have_triangle = 0;
    f_count++;
    return SM_TRIANGLE;
  }
  else if (have_finalized)
  {
    final_idx = t_idx[finalized_vertices[next_finalized]];
    have_finalized--; next_finalized++;
    return SM_FINALIZED;
  }
  return SM_ERROR;
}

SMreadWithSecondDerivatives::SMreadWithSecondDerivatives()
{
  // init of SMreadWithSecondDerivatives
  v_nor_f[0] = v_nor_f[1] = v_nor_f[2] = 0;

  t_pos_f[0] = t_pos_f[1] = t_pos_f[2] = 0;
  t_nor_f[0] = t_nor_f[1] = t_nor_f[2] = 0;

  smreadwithvertexnormals = 0;

  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;
}

SMreadWithSecondDerivatives::~SMreadWithSecondDerivatives()
{
  // clean-up for SMreadWithSecondDerivatives
  bb_min_f = 0;
  bb_max_f = 0;
  if (smreadwithvertexnormals) delete smreadwithvertexnormals;
}
