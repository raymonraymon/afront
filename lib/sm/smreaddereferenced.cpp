/*
===============================================================================

  FILE:  smreaddereferenced.cpp
  
  CONTENTS:
  
    see corresponding header file

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2007  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see corresponding header file
  
===============================================================================
*/
#include "smreaddereferenced.h"

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
  int index;
  short list_size;
  short list_alloc;
  SMtriangle** list;
} SMvertex;

typedef struct SMtriangle
{
  SMtriangle* buffer_next;  // used for efficient memory management
  SMvertex* vertices[3];
  char arrived_vertices;
} SMtriangle;

typedef hash_map<int, SMvertex*> my_vertex_hash;
typedef deque<SMtriangle*> my_triangle_queue;

my_vertex_hash* vertex_hash = 0;
my_triangle_queue* triangle_queue = 0;

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
  vertex->list_size = 0;
  vertex->buffer_next = 0;

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
  triangle->arrived_vertices = 0;

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

bool SMreadDereferenced::open(SMreader* smreader)
{
  if (smreader == 0)
  {
    return false;
  }
  this->smreader = smreader;

  nverts = smreader->nverts;
  nfaces = smreader->nfaces;

  v_count = 0;
  f_count = 0;

  bb_min_f = smreader->bb_min_f;
  bb_max_f = smreader->bb_max_f;

  vertex_hash = new my_vertex_hash;
  triangle_queue = new my_triangle_queue;

  initVertexBuffer(1024);
  initTriangleBuffer(2048);

  return true;
}

void SMreadDereferenced::close(bool close_file)
{
  // close of SMreadDereferenced
  smreader->close(close_file);

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

int SMreadDereferenced::read_input()
{
  int i;
  SMvertex* vertex;
  SMtriangle* triangle;
  my_vertex_hash::iterator hash_elements[3];

  SMevent event = smreader->read_element();
    
  if (event == SM_TRIANGLE)
  {
    // create triangle
    triangle = allocTriangle();
    // lookup the triangle's vertices
    for (i = 0; i < 3; i++)
    {
      // look for vertex in the hash
      hash_elements[i] = vertex_hash->find(smreader->t_idx[i]);
      // did we find the vertex in the hash  
      if (hash_elements[i] == vertex_hash->end())
      {
        // vertex has not arrived yet. so we create the vertex ...
        vertex = allocVertex();
        // ... insert it into hash ...
        vertex_hash->insert(my_vertex_hash::value_type(smreader->t_idx[i], vertex));
        // ... and mark it as not having arrived
        vertex->index = -2;
      }
      else
      {
        // use vertex found in hash
        vertex = (*hash_elements[i]).second;
      }
      // add vertex to triangle
      triangle->vertices[i] = vertex;
      // has this vertex arrived already
      if (vertex->index == -2)
      {
        // no ... we store the triangle in a list so we can tell it about the arrival later
        if (vertex->list)
        {
          if (vertex->list_size == vertex->list_alloc)
          {
            vertex->list_alloc += 4;
            vertex->list = (SMtriangle**)realloc(vertex->list,vertex->list_alloc*sizeof(SMtriangle*));
          }
        }
        else
        {
          vertex->list_alloc = 10;
          vertex->list = (SMtriangle**)malloc(vertex->list_alloc*sizeof(SMtriangle*));
          assert(vertex->list_size == 0);
        }
        vertex->list[vertex->list_size] = triangle;
      }
      else
      {
        // yes ... we can tell the triangle the triangle now
        triangle->arrived_vertices++;
      }
      // vertex has one more triangle
      vertex->list_size++;
    }
    // have all vertices of this triangle already arrived
    if (triangle->arrived_vertices == 3)
    {
      // then it is ready to be output
      triangle_queue->push_back(triangle);
    }
    // check for finalization of the triangle's vertices
    for (i = 0; i < 3; i++)
    {
      // look for vertex in the hash
      if (smreader->t_final[i])
      {
        // mark vertex as finalized
        triangle->vertices[i]->buffer_next = triangle->vertices[i];
        // remove vertex from hash
        vertex_hash->erase(hash_elements[i]);
      }
    }
  }
  else if (event == SM_VERTEX)
  {
    // look for vertex in the vertex hash
    hash_elements[0] = vertex_hash->find(smreader->v_idx);
    // did we find the vertex in the hash  
    if (hash_elements[0] == vertex_hash->end())
    {
      // vertex had not been referenced yet. so we create the vertex ...
      vertex = allocVertex();
      // ... insert it into hash
      vertex_hash->insert(my_vertex_hash::value_type(smreader->v_idx, vertex));
    }
    else
    {
      // use vertex found in hash
      vertex = vertex = (*hash_elements[0]).second;
      // make sure it is not yet marked as arrived
      assert(vertex->index == -2);
      // make sure it has at least one triangle in its list
      assert(vertex->list_size > 0);
    }
    // copy coordinates
    VecCopy3fv(vertex->v, smreader->v_pos_f);
    // and mark the vertex as having arrived
    vertex->index = -1;
    // did we already have triangles
    if (vertex->list_size)
    {
      // then loop over them to tell them about the arrival
      for (i = 0; i < vertex->list_size; i++)
      {
        // get the triangle
        triangle = vertex->list[i];
        // increment the arrival counter
        triangle->arrived_vertices++;
        // is this the third vertex to arrive
        if (triangle->arrived_vertices == 3)
        {
          // then we can output it
          triangle_queue->push_back(triangle);
        }
      }
    }
    // in post order the arrival of a vertex marks its finalization
    if (smreader->post_order)
    {
      // mark vertex as finalized
      vertex->buffer_next = vertex;
      // so we remove it from hash
      vertex_hash->erase(hash_elements[0]);
    }
  }
  else if (event == SM_FINALIZED)
  {
    // look for vertex in the vertex hash
    hash_elements[0] = vertex_hash->find(smreader->v_idx);
    // did we find the vertex in the hash  
    if (hash_elements[0] == vertex_hash->end())
    {
      // vertex was not found in hash
      fprintf(stderr, "FATAL ERROR: finalized vertex %d not in hash. corrupt streaming mesh.\n", smreader->v_idx);
      return -1;
    }
    else
    {
      // use vertex found in hash
      vertex = vertex = (*hash_elements[0]).second;
    }
    // mark vertex as finalized
    vertex->buffer_next = vertex;
    // so we remove it from hash
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
      // ... mark this vertex as finalized ...
      vertex->buffer_next = vertex;
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

int SMreadDereferenced::read_triangle()
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

  assert(triangle->arrived_vertices == 3);
  
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
    // is this vertex to be finalized 
    if (vertex->list_size == 1 && vertex->buffer_next) // if this is the last triangle and the vertex was finalized
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
  }

  deallocTriangle(triangle);

  return 1;
}

SMevent SMreadDereferenced::read_element()
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

SMevent SMreadDereferenced::read_event()
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

SMreadDereferenced::SMreadDereferenced()
{
  // init of SMreadDereferenced
  t_pos_f[0] = t_pos_f[1] = t_pos_f[2] = 0;

  smreader = 0;

  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;
}

SMreadDereferenced::~SMreadDereferenced()
{
  // clean-up for SMreadDereferenced
  bb_min_f = 0;
  bb_max_f = 0;
  if (smreader) delete smreader;
}
