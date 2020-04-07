/*
===============================================================================

  FILE:  smreadbuffered.cpp
  
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
#include "smreadbuffered.h"

#include <stdio.h>

#include "vec3fv.h"
#include "vec3iv.h"

#include "hash_map.h"

#include "dynamicqueue.h"

#define PRINT_CONTROL_OUTPUT
#undef PRINT_CONTROL_OUTPUT

// data structures used for streaming re-ordering 

struct SMtriangle;

typedef struct SMvertex
{
  SMvertex* buffer_next; // used for efficient memory management
  int index;
  bool finalized;
  float v[3];
  // triangles
  int incoming_size;
  int incoming_alloc;
  SMtriangle** incoming;
} SMvertex;

typedef struct SMtriangle
{
  SMtriangle* buffer_next; // used for efficient memory management (and by the waiting queue)
  int dirty;               // to dirty triangle in the buffer (so they are only checked once)
  SMvertex* vertices[5];
} SMtriangle;

typedef hash_map<int, SMvertex*> my_hash;

static my_hash* vertex_hash; // for matching vertices

static DynamicQueue* waiting_queue; // for triangles ready for output

static SMvertex* cache[6]; // for simulating the cache (and avoid '(i+1)%3'-style computations)

static SMvertex** vertices;

// statistics

#ifdef PRINT_CONTROL_OUTPUT
static int triangle_buffer_maxsize = 0;
static int vertex_buffer_maxsize = 0;
#endif

// efficient memory allocation for vertices

static int vertex_buffer_size = 0;
static int vertex_buffer_alloc = 1024;
static SMvertex* vertex_buffer_next = 0;

static int initVertexBuffer(int size)
{
  vertex_buffer_next = (SMvertex*)malloc(sizeof(SMvertex)*size);

  if (vertex_buffer_next == 0)
  {
    fprintf(stderr,"malloc for vertex buffer failed\n");
    return 0;
  }
  for (int i = 0; i < size; i++)
  {
    vertex_buffer_next[i].buffer_next = &(vertex_buffer_next[i+1]);
    vertex_buffer_next[i].incoming = 0;
  }
  vertex_buffer_next[size-1].buffer_next = 0;
  vertex_buffer_alloc = size;
  vertex_buffer_size = 0;
  return 1;
}

static SMvertex* allocVertex()
{
  if (vertex_buffer_next == 0)
  {
    vertex_buffer_next = (SMvertex*)malloc(sizeof(SMvertex)*vertex_buffer_alloc);
    if (vertex_buffer_next == 0)
    {
      fprintf(stderr,"malloc for vertex buffer failed\n");
      return 0;
    }
    for (int i = 0; i < vertex_buffer_alloc; i++)
    {
      vertex_buffer_next[i].buffer_next = &(vertex_buffer_next[i+1]);
      vertex_buffer_next[i].incoming = 0;
    }
    vertex_buffer_next[vertex_buffer_alloc-1].buffer_next = 0;
    vertex_buffer_alloc = 2*vertex_buffer_alloc;
  }
  // get pointer to next available vertex
  SMvertex* vertex = vertex_buffer_next;
  vertex_buffer_next = vertex->buffer_next;

  vertex->finalized = false;
  vertex->index = -1;
  if (vertex->incoming == 0)
  {
    vertex->incoming = (SMtriangle**)malloc(sizeof(SMtriangle*)*10);
    vertex->incoming_alloc = 10;
  }
  vertex->incoming_size = 0;

  vertex_buffer_size++;
  
#ifdef PRINT_CONTROL_OUTPUT
  if (vertex_buffer_size > vertex_buffer_maxsize) vertex_buffer_maxsize = vertex_buffer_size;
#endif

  return vertex;
}

static void deallocVertex(SMvertex* vertex)
{
  vertex->buffer_next = vertex_buffer_next;
  vertex_buffer_next = vertex;
  vertex_buffer_size--;
}

// efficient memory allocation for triangles

static int triangle_buffer_size = 0;
static int triangle_buffer_alloc = 1024;
static SMtriangle* triangle_buffer_next = 0;

static int initTriangleBuffer(int alloc)
{
  triangle_buffer_next = (SMtriangle*)malloc(sizeof(SMtriangle)*alloc);

  if (triangle_buffer_next == 0)
  {
    fprintf(stderr,"malloc for triangle buffer failed\n");
    return 0;
  }
  for (int i = 0; i < alloc; i++)
  {
    triangle_buffer_next[i].buffer_next = &(triangle_buffer_next[i+1]);
  }
  triangle_buffer_next[alloc-1].buffer_next = 0;
  triangle_buffer_alloc = alloc;
  triangle_buffer_size = 0;
  return 1;
}

static SMtriangle* allocTriangle()
{
  if (triangle_buffer_next == 0)
  {
    triangle_buffer_next = (SMtriangle*)malloc(sizeof(SMtriangle)*triangle_buffer_alloc);
    if (triangle_buffer_next == 0)
    {
      fprintf(stderr,"malloc for triangle buffer failed\n");
      return 0;
    }
    for (int i = 0; i < triangle_buffer_alloc; i++)
    {
      triangle_buffer_next[i].buffer_next = &(triangle_buffer_next[i+1]);
    }
    triangle_buffer_next[triangle_buffer_alloc-1].buffer_next = 0;
    triangle_buffer_alloc = 2*triangle_buffer_alloc;
  }
  // get pointer to next available vertex
  SMtriangle* triangle = triangle_buffer_next;
  triangle_buffer_next = triangle->buffer_next;

  triangle->dirty = -1;
  triangle_buffer_size++;
  
#ifdef PRINT_CONTROL_OUTPUT
  if (triangle_buffer_size > triangle_buffer_maxsize) triangle_buffer_maxsize = triangle_buffer_size;
#endif

  return triangle;
}

static void deallocTriangle(SMtriangle* triangle)
{
  triangle->buffer_next = triangle_buffer_next;
  triangle_buffer_next = triangle;
  triangle_buffer_size--;
}

static void addToVertex(SMtriangle* triangle, SMvertex* vertex)
{
  if (vertex->incoming_size == vertex->incoming_alloc)
  {
    vertex->incoming_alloc += 2;
    vertex->incoming = (SMtriangle**)realloc(vertex->incoming, sizeof(SMtriangle*)*vertex->incoming_alloc);
  }
  vertex->incoming[vertex->incoming_size] = triangle;
  vertex->incoming_size++;
}

static void addToVertices(SMtriangle* triangle)
{
  addToVertex(triangle, triangle->vertices[0]);
  addToVertex(triangle, triangle->vertices[1]);
  addToVertex(triangle, triangle->vertices[2]);
}

static void removeFromVertex(SMtriangle* triangle, SMvertex* vertex)
{
  for (int i = 0; i < vertex->incoming_size; i++)
  {
    if (vertex->incoming[i] == triangle)
    {
      vertex->incoming_size--;
      vertex->incoming[i] = vertex->incoming[vertex->incoming_size];
      return;
    }
  }
  fprintf(stderr,"ERROR: did not find triangle in incoming\n");
  exit(0);
}

static void removeFromVertices(SMtriangle* triangle)
{
  removeFromVertex(triangle, triangle->vertices[0]);
  removeFromVertex(triangle, triangle->vertices[1]);
  removeFromVertex(triangle, triangle->vertices[2]);
}

bool SMreadBuffered::open(SMreader* smreader, int max_delay)
{
  if (smreader == 0 )
  {
    return false;
  }
  this->smreader = smreader;
  this->max_delay = max_delay;

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"maximal triangle delay is %d\n",max_delay);
#endif
  
  nverts = smreader->nverts;
  nfaces = smreader->nfaces;

  v_count = 0;
  f_count = 0;

  bb_min_f = smreader->bb_min_f;
  bb_max_f = smreader->bb_max_f;

  vertex_hash = new my_hash;
  waiting_queue = new DynamicQueue();

  initVertexBuffer(1024);
  initTriangleBuffer(2048);

  cache[0] = 0;
  cache[1] = 0;
  cache[2] = 0;
  cache[3] = 0;
  cache[4] = 0;
  cache[5] = 0;

  return true;
}

void SMreadBuffered::close(bool close_file)
{
  // close of SMreadBuffered
  smreader->close(close_file);

  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;

  delete vertex_hash;
  delete waiting_queue;

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

static SMtriangle* find_triangle(int dirty)
{
  int i,j,k;
  SMtriangle* triangle;
  SMvertex** verts;
  int ccx_age = 2000000000;
  SMtriangle* ccx = 0; // triangle with two cache vertices
  int cax_age = -1;
  SMtriangle* caf = 0; // triangle with one cache vertex and two active vertices of which one will be finalized
  SMtriangle* cax = 0; // triangle with one cache vertex and one active vertex
  SMtriangle* caa = 0; // triangle with one cache vertex and two active vertices

  // check if there are good triangles around cache entry i

  for (i = 0; i < 3; i++)
  {
    if (cache[i] && cache[i]->incoming_size) // does this cache vertex has triangles
    {
      if (cache[i]->incoming_size == 1 && cache[i]->finalized) // does the next triangle finalize this vertex
      {
        return cache[i]->incoming[0]; // do it. finalized vertices can not cause future cache misses
      }
      for (j = 0; j < cache[i]->incoming_size; j++) // iterate over the triangles of cache vertex i
      {
        triangle = cache[i]->incoming[j];
        if (triangle->dirty != dirty) // if not checked yet
        {
          triangle->dirty = dirty; // mark as checked
          verts = triangle->vertices; // get its vertices
          for (k = 0; k < 3; k++) // which of this triangle's vertices is the cache vertex i
          {
            if (cache[i] == verts[k]) // is it triangle vertex k
            {
              if (cache[i+1] == verts[k+2]) //is triangle vertex k+2 also in cache
              {
                if (verts[k+1]->incoming_size == 1 && verts[k+1]->finalized)
                {
                  return verts[k+1]->incoming[0]; // use it. it finalizes the third triangle vertex k+1
                }
                if (verts[k]->index < verts[k+2]->index)
                {
                  if (verts[k]->index < ccx_age) // continue with the oldest
                  {
                    ccx_age = verts[k]->index;
                    ccx = triangle;
                  }
                }
                else
                {
                  if (verts[k+2]->index < ccx_age) // continue with the oldest
                  {
                    ccx_age = verts[k+2]->index;
                    ccx = triangle;
                  }
                }
              }
              else if (cache[i+2] == verts[k+1]) // else if v1 of triangle is also in cache
              {
                if (verts[k+2]->incoming_size == 1 && verts[k+2]->finalized)
                {
                  return verts[k+2]->incoming[0];
                }
                if (verts[k]->index < verts[k+1]->index)
                {
                  if (verts[k]->index < ccx_age) // continue with the oldest
                  {
                    ccx_age = verts[k]->index;
                    ccx = triangle;
                  }
                }
                else
                {
                  if (verts[k+1]->index < ccx_age) // continue with the oldest
                  {
                    ccx_age = verts[k+1]->index;
                    ccx = triangle;
                  }
                }
              }
              else if (verts[k+2]->index >= 0) // else if v2 of triangle is already active
              {
                if (verts[k+1]->index >= 0) // if v1 is also active
                {
                  if (verts[k+1]->incoming_size == 1 && verts[k+1]->finalized) // if this also finalizes v1
                  {
                    caf = triangle;
                  }
                  else if (verts[k+2]->incoming_size == 1 && verts[k+2]->finalized) // if this also finalizes v2
                  {
                    caf = triangle;
                  }
                  else
                  {
                    caa = triangle;
                  }
                }
                else
                {
                  if (verts[k+2]->index > cax_age) // continue with the most recent
                  {
                    cax_age = verts[k+2]->index;
                    cax = triangle;
                  }
                }
              }
              else if (verts[k+1]->index >= 0) // else if v1 of triangle is already active
              {
                if (verts[k+1]->index > cax_age) // continue with the most recent
                {
                  cax_age = verts[k+1]->index;
                  cax = triangle;
                }
              }
              break;
            }
          }
        }
      }
    }
  }

  if (ccx) return ccx; // first double cache hits
  else if (caf) return caf; // then cache hits that finalize a vertex
  else if (cax) return cax; // then cache hits with only one other vertex active (to avoid bias towards joins)
  else if (caa) return caa; // only then cache hits with both other vertices active 
  else if (cache[0] && cache[0]->incoming_size) return cache[0]->incoming[0];
  else if (cache[1] && cache[1]->incoming_size) return cache[1]->incoming[0];
  else if (cache[2] && cache[2]->incoming_size) return cache[2]->incoming[0];
  else return 0;
}

void SMreadBuffered::fill_triangle_buffer()
{
  SMevent event;

  while (waiting_queue->size() < max_delay)
  {
    event = smreader->read_element();
    if (event == SM_TRIANGLE)
    {
      int i;
      my_hash::iterator hash_elements[3];

      SMtriangle* triangle = allocTriangle();

      // get vertices from hash
      for (i = 0; i < 3; i++)
      {
        hash_elements[i] = vertex_hash->find(smreader->t_idx[i]);
        if (hash_elements[i] == vertex_hash->end())
        {
          fprintf(stderr,"FATAL ERROR: vertex not in hash. need pre-order mesh\n");
          exit(0);
        }
        else
        {
          triangle->vertices[i] = (*hash_elements[i]).second;
        }
      }

      triangle->vertices[3] = triangle->vertices[0];
      triangle->vertices[4] = triangle->vertices[1];
  
      addToVertices(triangle);

      // check for finalization
      for (i = 0; i < 3; i++)
      {
        if (smreader->t_final[i])
        {
          triangle->vertices[i]->finalized = true;
          vertex_hash->erase(hash_elements[i]);
        }
      }

      // add triangle to buffer
      waiting_queue->addElement(triangle);
    }
    else if (event == SM_VERTEX)
    {
      SMvertex* vertex = allocVertex();
      VecCopy3fv(vertex->v, smreader->v_pos_f);
      vertex_hash->insert(my_hash::value_type(smreader->v_idx, vertex));
    }
    else if (event == SM_FINALIZED)
    {
      my_hash::iterator hash_element = vertex_hash->find(smreader->final_idx);
      if (hash_element == vertex_hash->end())
      {
        fprintf(stderr,"FATAL ERROR: vertex not in hash. need pre-order mesh\n");
        exit(0);
      }
      else
      {
        ((*hash_element).second)->finalized = true;
        vertex_hash->erase(hash_element);
      }
    }
    else
    {
      max_delay = 0;
    }
  }
}

int SMreadBuffered::read_triangle_delayed()
{
  fill_triangle_buffer();

  SMtriangle* triangle = find_triangle(f_count);

  if (triangle)
  {
    waiting_queue->removeElement(triangle);
  }
  else
  {
    if (waiting_queue->elements())
    {
      triangle = (SMtriangle*)waiting_queue->getAndRemoveFirstElement();
    }
    else
    {
      return 0;
    }
  }

  removeFromVertices(triangle);

  vertices = triangle->vertices;

  for (int i = 0; i < 3; i++)
  {
    if (vertices[i]->index == -1)
    {
      vertices[i]->index = v_count + have_new;
      new_vertices[have_new] = i;
      have_new++;
    }
    t_idx[i] = vertices[i]->index;
    if (vertices[i]->finalized && vertices[i]->incoming_size == 0)
    {
      t_final[i] = true;
      deallocVertex(vertices[i]);
      cache[i] = 0;
      cache[i+3] = 0;
      finalized_vertices[have_finalized] = i;
      have_finalized++;
    }
    else
    {
      t_final[i] = false;
      cache[i] = vertices[i];
      cache[i+3] = vertices[i];;
    }
  }
  deallocTriangle(triangle);
  have_triangle = 1;

  return 1;
}

SMevent SMreadBuffered::read_element()
{
  if (have_triangle == 0)
  {
    next_new = 0;
    have_finalized = next_finalized = 0;
    if (read_triangle_delayed() == 0)
    {
      if (nverts != -1 && v_count != nverts)
      {
        fprintf(stderr,"WARNING: wrong vertex count: v_count (%d) != nverts (%d)\n", v_count, nverts);
      }
      nverts = v_count;
      if (nfaces != -1 && f_count != nfaces)
      {
        fprintf(stderr,"WARNING: wrong face count: f_count (%d) != nfaces (%d)\n", f_count, nfaces);
      }
      nfaces = f_count;
      have_triangle = -1;
      return SM_EOF;
    }
  }
  if (have_new)
  {
    v_idx = t_idx[new_vertices[next_new]];
    VecCopy3fv(v_pos_f, vertices[new_vertices[next_new]]->v);
    have_new--; next_new++;
    v_count++;
    return SM_VERTEX;
  }
  else if (have_triangle > 0)
  {
    have_triangle = 0;
    f_count++;
    return SM_TRIANGLE;
  }
  else if (have_triangle == -1)
  {
    return SM_EOF;
  }
  return SM_ERROR;
}

SMevent SMreadBuffered::read_event()
{
  if (have_finalized)
  {
    final_idx = finalized_vertices[next_finalized];
    have_finalized--; next_finalized++;
    return SM_FINALIZED;
  }
  else
  {
    return read_element();
  }
}

SMreadBuffered::SMreadBuffered()
{
  // init of SMreadBuffered
  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;
  smreader = 0;
  max_delay = -1;
}

SMreadBuffered::~SMreadBuffered()
{
  // clean-up for SMreadBuffered
  bb_min_f = 0;
  bb_max_f = 0;
  if (smreader) delete smreader;
}
