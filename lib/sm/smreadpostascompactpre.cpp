/*
===============================================================================

  FILE:  smreadpostascompactpre.cpp
  
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
#include "smreadpostascompactpre.h"

#include <stdio.h>

#include "vec3fv.h"
#include "vec3iv.h"

#include "hash_map.h"

#include "dynamicvector.h"

struct SMtriangle;

typedef struct SMvertex
{
  SMvertex* buffer_next;    // used for efficient memory management
  float v[3];
  int index;
  int list_size;
  int list_alloc;
  SMtriangle** list;
} SMvertex;

typedef struct SMtriangle
{
  int dynamicvector;        // used by the dynamicvector
  SMtriangle* buffer_next;  // used for efficient memory management
  SMvertex* vertices[3];
  int ready;
} SMtriangle;

typedef hash_map<int, SMvertex*> my_hash;

static my_hash* vertex_hash;

static DynamicVector* waiting_area;
static DynamicVector* output_triangles;

static int vertex_buffer_size;
static int triangle_buffer_size;

static int vertex_buffer_maxsize;
static int triangle_buffer_maxsize;

// efficient memory allocation

static SMvertex* vertex_buffer_next;
static int vertex_buffer_alloc;
static int vertex_buffer_increment;
static int vertex_buffer_limit;
static int vertex_buffer_number_chunks;
static SMvertex* vertex_buffer_chunks[16];

static void initVertexBuffer(int alloc, int limit)
{
  vertex_buffer_next = 0;
  vertex_buffer_alloc = 0;
  vertex_buffer_increment = alloc;
  vertex_buffer_limit = limit;
  vertex_buffer_number_chunks = 0;
}

static SMvertex* allocVertex()
{
  if (vertex_buffer_next == 0)
  {
    int increment;
    if (vertex_buffer_limit && (vertex_buffer_alloc + vertex_buffer_increment >= vertex_buffer_limit))
    {
      if (vertex_buffer_alloc == vertex_buffer_limit)
      {
        fprintf(stderr,"FATAL ERROR: reached vertex buffer limit\n");
        return 0;
      }
      increment = vertex_buffer_limit - vertex_buffer_alloc;
    }
    else
    {
      increment = vertex_buffer_increment;
      vertex_buffer_increment*=2;
    }
    vertex_buffer_next = (SMvertex*)malloc(sizeof(SMvertex)*increment);
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
    vertex_buffer_alloc += increment;
    for (int i = 0; i < increment; i++)
    {
      vertex_buffer_next[i].buffer_next = &(vertex_buffer_next[i+1]);
      vertex_buffer_next[i].list = 0;
    }
    vertex_buffer_next[increment-1].buffer_next = 0;
  }
  // get pointer to next available vertex
  SMvertex* vertex = vertex_buffer_next;
  vertex_buffer_next = vertex->buffer_next;

  // initialize data fields
  vertex->index = -1;
  if (vertex->list == 0)
  {
    vertex->list = (SMtriangle**)malloc(sizeof(SMtriangle*)*10);
    vertex->list_alloc = 10;
  }
  vertex->list_size = 0;

  // keep track of maxsize
  vertex_buffer_size++; if (vertex_buffer_size > vertex_buffer_maxsize) vertex_buffer_maxsize = vertex_buffer_size;

  return vertex;
}

static void deallocVertex(SMvertex* vertex)
{
  vertex->buffer_next = vertex_buffer_next;
  vertex_buffer_next = vertex;
  vertex_buffer_size--;
}

static void destroyVertexBuffer()
{
  for (int i = 0; i < vertex_buffer_number_chunks; i++)
  {
    free(vertex_buffer_chunks[i]);
  }
  vertex_buffer_number_chunks = 0;
  vertex_buffer_next = 0;
  vertex_buffer_alloc = 0;
}

static SMtriangle* triangle_buffer_next;
static int triangle_buffer_alloc;
static int triangle_buffer_increment;
static int triangle_buffer_limit;
static int triangle_buffer_number_chunks;
static SMtriangle* triangle_buffer_chunks[16];

static void initTriangleBuffer(int alloc, int limit)
{
  triangle_buffer_next = 0;
  triangle_buffer_alloc = 0;
  triangle_buffer_increment = alloc;
  triangle_buffer_limit = limit;
  triangle_buffer_number_chunks = 0;
}

static SMtriangle* allocTriangle()
{
  if (triangle_buffer_next == 0)
  {
    int increment;
    if (triangle_buffer_limit && (triangle_buffer_alloc + triangle_buffer_increment >= triangle_buffer_limit))
    {
      if (triangle_buffer_alloc == triangle_buffer_limit)
      {
        fprintf(stderr,"FATAL ERROR: reached triangle buffer limit\n");
        return 0;
      }
      increment = triangle_buffer_limit - triangle_buffer_alloc;
    }
    else
    {
      increment = triangle_buffer_increment;
      triangle_buffer_increment*=2;
    }
    triangle_buffer_next = (SMtriangle*)malloc(sizeof(SMtriangle)*increment);
    if (triangle_buffer_next == 0)
    {
      fprintf(stderr,"FATAL ERROR: malloc for triangle buffer failed\n");
      return 0;
    }
    for (int i = 0; i < increment; i++)
    {
      triangle_buffer_next[i].buffer_next = &(triangle_buffer_next[i+1]);
    }
    triangle_buffer_next[increment-1].buffer_next = 0;
  }
  // get pointer to next available triangle
  SMtriangle* triangle = triangle_buffer_next;
  triangle_buffer_next = triangle->buffer_next;
 
  // initialize data fields
  triangle->ready = 0;

  // keep track of maxsize
  triangle_buffer_size++; if (triangle_buffer_size > triangle_buffer_maxsize) triangle_buffer_maxsize = triangle_buffer_size;

  return triangle;
}

static void deallocTriangle(SMtriangle* triangle)
{
  triangle->buffer_next = triangle_buffer_next;
  triangle_buffer_next = triangle;
  triangle_buffer_size--;
}

static void destroyTriangleBuffer()
{
  for (int i = 0; i < triangle_buffer_number_chunks; i++)
  {
    free(triangle_buffer_chunks[i]);
  }
  triangle_buffer_number_chunks = 0;
  triangle_buffer_next = 0;
  triangle_buffer_alloc = 0;
}

bool SMreadPostAsCompactPre::open(SMreader* smreader, bool preserve_order, int limit_buffer_size)
{
  if (smreader == 0 || !smreader->post_order)
  {
    return false;
  }
  this->smreader = smreader;
  this->preserve_order = preserve_order;

  nverts = smreader->nverts;
  nfaces = smreader->nfaces;

  v_count = 0;
  f_count = 0;

  bb_min_f = smreader->bb_min_f;
  bb_max_f = smreader->bb_max_f;

  vertex_hash = new my_hash;

  if (preserve_order)
  {
    waiting_area = new DynamicVector();
  }
  else
  {
    waiting_area = 0;
  }
  output_triangles = new DynamicVector();

  initVertexBuffer(1024, 0);
  initTriangleBuffer(2048, limit_buffer_size);

  return true;
}

void SMreadPostAsCompactPre::close(bool close_file)
{
  // close of SMreadPostAsCompactPre
  smreader->close(close_file);

  delete vertex_hash;

  if (preserve_order) delete waiting_area;
  delete output_triangles;

  destroyVertexBuffer();
  destroyTriangleBuffer();

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

int SMreadPostAsCompactPre::read_postorder_input()
{
  int i;
  SMvertex* vertex;
  SMtriangle* triangle;
  my_hash::iterator hash_element;

  SMevent event = smreader->read_element();
    
  if (event == SM_TRIANGLE)
  {
    // create triangle
    triangle = allocTriangle();
    // insert triangle into waiting area
    if (preserve_order) waiting_area->addElement(triangle);
    // process the triangle's vertices
    for (i = 0; i < 3; i++)
    {
      // look for vertex in the hash
      hash_element = vertex_hash->find(smreader->t_idx[i]);
      // did we find the vertex in the hash  
      if (hash_element == vertex_hash->end())
      {
        // create vertex
        vertex = allocVertex();
        // insert vertex into hash
        vertex_hash->insert(my_hash::value_type(smreader->t_idx[i], vertex));
      }
      else
      {
        // use vertex found in hash
        vertex = (*hash_element).second;
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
  }
  else if (event == SM_VERTEX)
  {
    // look for vertex in the vertex hash
    hash_element = vertex_hash->find(smreader->v_idx);
    // we inform the user if a vertex was not in the hash
    if (hash_element == vertex_hash->end())
    {
      fprintf(stderr, "WARNING: post-order vertex not used by any triangle. skipping ...\n");
      return 1;
    }
    // use vertex found in hash
    vertex = vertex = (*hash_element).second;
    // copy coordinates
    VecCopy3fv(vertex->v, smreader->v_pos_f);
    // finalizing vertex may make some triangles eligible for output
    for (i = 0; i < vertex->list_size; i++)
    {
      triangle = vertex->list[i];
      triangle->ready++;
      if (triangle->ready == 3)
      {
        if (preserve_order) // we have to preserve the triangle order
        {
          if (triangle == waiting_area->getFirstElement())
          {
            waiting_area->removeFirstElement();
            output_triangles->addElement(triangle);

            // behind this triangle additional ready triangles may be waiting

            while (waiting_area->size())
            {
              triangle = (SMtriangle*)waiting_area->getFirstElement();
              if (triangle->ready == 3)
              {
                waiting_area->removeFirstElement();
                output_triangles->addElement(triangle);
              }
              else
              {
                break;
              }
            }
          }
        }
        else // we do not have to preserve the triangle order
        {
          output_triangles->addElement(triangle);
        }
      }
    }
    // remove vertex from hash
    vertex_hash->erase(hash_element);
  }
  else if (event == SM_EOF)
  {
    if (vertex_hash->size())
    {
      fprintf(stderr, "FATAL ERROR: some outstanding vertices when reaching SM_EOF\n");
      return -1;
    }
    else
    {
      return 0;
    }
  }
  return 1;
}

int SMreadPostAsCompactPre::read_triangle()
{
  int i;
  SMvertex* vertex;
  SMtriangle* triangle;

  // read from postorder input until we have a compact preorder triangle for output (or fail)

  while (!output_triangles->size())
  {
    int ok = read_postorder_input();
    if (ok <= 0)
    {
      return ok;
    }
  }

  // get next output triangle

  triangle = (SMtriangle*)output_triangles->getAndRemoveFirstElement();
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
    // is this vertex to be finalized?
    if (vertex->list_size == 1)
    {
      t_final[i] = true;
      finalized_vertices[have_finalized++] = i;
      deallocVertex(vertex); // it can still be used until the next vertex is alloced.
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
  }

  deallocTriangle(triangle);

  return 1;
}

SMevent SMreadPostAsCompactPre::read_element()
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
    VecCopy3fv(v_pos_f,t_pos_f[new_vertices[next_new]]);
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

SMevent SMreadPostAsCompactPre::read_event()
{
  if ((have_new + have_triangle + have_finalized) == 0)
  {
    return read_element();
  }
  if (have_new)
  {
    v_idx = t_idx[new_vertices[next_new]];
    VecCopy3fv(v_pos_f, t_pos_f[new_vertices[next_new]]);
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

SMreadPostAsCompactPre::SMreadPostAsCompactPre()
{
  // init of SMreadPostAsCompactPre
  t_pos_f[0] = t_pos_f[1] = t_pos_f[2] = 0;

  smreader = 0;
  preserve_order = true;

  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;
}

SMreadPostAsCompactPre::~SMreadPostAsCompactPre()
{
  // clean-up for SMreadPostAsCompactPre
  bb_min_f = 0;
  bb_max_f = 0;
  if (smreader) delete smreader;
}
