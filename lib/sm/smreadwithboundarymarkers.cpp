/*
===============================================================================

  FILE:  smreadwithboundarymarkers.cpp
  
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
#include "smreadwithboundarymarkers.h"

#include <stdlib.h>
#include <stdio.h>

#include "hash_map.h"
#include "deque.h"

#include "vec3fv.h"

#define PRINT_CONTROL_OUTPUT

struct SMtriangle;

typedef struct SMvertex
{
  SMvertex* buffer_next;  // used for efficient memory management
  float pos_f[3];
  int index;
  bool boundary;
  int list_alloc;
  int list_size;
  SMtriangle** list;
} SMvertex;

typedef struct SMtriangle
{
  SMtriangle* buffer_next;  // used for efficient memory management
  SMvertex* vertices[3];
  int ready;
} SMtriangle;

typedef hash_map<int, SMvertex*> my_hash;
typedef deque<SMtriangle*> my_queue;

static my_hash* vertex_hash;
static my_queue* triangle_queue;

// efficient memory allocation

static SMvertex* vertex_buffer_next;
static int vertex_buffer_alloc;
static int vertex_buffer_increment;
static int vertex_buffer_initial_increment;
static int vertex_buffer_number_chunks;
static SMvertex* vertex_buffer_chunks[16];

#ifdef PRINT_CONTROL_OUTPUT
static int vertex_buffer_size;
static int vertex_buffer_maxsize;
#endif

static void initVertexBuffer(int alloc)
{
  vertex_buffer_next = 0;
  vertex_buffer_alloc = 0;
  vertex_buffer_increment = alloc;
  vertex_buffer_initial_increment = alloc;
  vertex_buffer_number_chunks = 0;
#ifdef PRINT_CONTROL_OUTPUT
  vertex_buffer_size = 0;
  vertex_buffer_maxsize = 0;
#endif
}

static SMvertex* allocVertex()
{
  if (vertex_buffer_next == 0)
  {
    vertex_buffer_next = (SMvertex*)malloc(sizeof(SMvertex)*vertex_buffer_increment);
    if (vertex_buffer_next == 0)
    {
      fprintf(stderr,"FATAL ERROR: malloc for vertex buffer with increment %d failed\n", vertex_buffer_increment);
      return 0;
    }
    else
    {
      vertex_buffer_chunks[vertex_buffer_number_chunks] = vertex_buffer_next;
      vertex_buffer_number_chunks++;
    }
    vertex_buffer_alloc += vertex_buffer_increment;
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
  vertex->index = -1;
  vertex->boundary = false;
  vertex->list_size = 0;

  if (vertex->list == 0)
  {
    vertex->list_alloc = 10;
    vertex->list = (SMtriangle**)malloc(sizeof(SMtriangle*)*10);
    if (vertex->list == 0)
    {
      fprintf(stderr,"FATAL ERROR: malloc for vertex->list with %d elements failed\n", 10);
      exit(1);
    }
  }

#ifdef PRINT_CONTROL_OUTPUT
  // keep track of maxsize
  vertex_buffer_size++; if (vertex_buffer_size > vertex_buffer_maxsize) vertex_buffer_maxsize = vertex_buffer_size;
#endif

  return vertex;
}

static void deallocVertex(SMvertex* vertex)
{
  vertex->buffer_next = vertex_buffer_next;
  vertex_buffer_next = vertex;
#ifdef PRINT_CONTROL_OUTPUT
  vertex_buffer_size--;
#endif
}

static void destroyVertexBuffer()
{
  for (int i = vertex_buffer_number_chunks-1; i >= 0; i--)
  {
    vertex_buffer_increment /= 2;
    for (int j = 0; j < vertex_buffer_increment; j++)
    {
      if (vertex_buffer_chunks[i][j].list) free (vertex_buffer_chunks[i][j].list);
    }
    free(vertex_buffer_chunks[i]);
  }
  vertex_buffer_number_chunks = 0;
  vertex_buffer_next = 0;
  vertex_buffer_alloc = 0;
#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr, "vertex_buffer_size %d vertex_buffer_maxsize %d\n", vertex_buffer_size, vertex_buffer_maxsize);
#endif
}

static SMtriangle* triangle_buffer_next;
static int triangle_buffer_alloc;
static int triangle_buffer_increment;
static int triangle_buffer_initial_increment;
static int triangle_buffer_number_chunks;
static SMtriangle* triangle_buffer_chunks[16];

#ifdef PRINT_CONTROL_OUTPUT
static int triangle_buffer_size;
static int triangle_buffer_maxsize;
#endif

static void initTriangleBuffer(int alloc)
{
  triangle_buffer_next = 0;
  triangle_buffer_alloc = 0;
  triangle_buffer_increment = alloc;
  triangle_buffer_initial_increment = alloc;
  triangle_buffer_number_chunks = 0;
#ifdef PRINT_CONTROL_OUTPUT
  triangle_buffer_size = 0;
  triangle_buffer_maxsize = 0;
#endif
}

static SMtriangle* allocTriangle()
{
  if (triangle_buffer_next == 0)
  {
    triangle_buffer_next = (SMtriangle*)malloc(sizeof(SMtriangle)*triangle_buffer_increment);
    if (triangle_buffer_next == 0)
    {
      fprintf(stderr,"FATAL ERROR: malloc for triangle buffer with increment %d failed\n", triangle_buffer_increment);
      exit(1);
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
  triangle->ready = 0;

#ifdef PRINT_CONTROL_OUTPUT
  // keep track of maxsize
  triangle_buffer_size++; if (triangle_buffer_size > triangle_buffer_maxsize) triangle_buffer_maxsize = triangle_buffer_size;
#endif

  return triangle;
}

static void deallocTriangle(SMtriangle* triangle)
{
  triangle->buffer_next = triangle_buffer_next;
  triangle_buffer_next = triangle;
#ifdef PRINT_CONTROL_OUTPUT
  triangle_buffer_size--;
#endif
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
#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr, "triangle_buffer_size %d triangle_buffer_maxsize %d\n", triangle_buffer_size, triangle_buffer_maxsize);
#endif
}

bool SMreadWithBoundaryMarkers::open(SMreader* smreader)
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

  vertex_hash = new my_hash;
  triangle_queue = new my_queue();

  initVertexBuffer(1024);
  initTriangleBuffer(2048);

  num_non_pre_order = 0;

  return true;
}

void SMreadWithBoundaryMarkers::close(bool close_file)
{
  smreader->close(close_file);

  v_count = -1;
  f_count = -1;

  delete vertex_hash;
  delete triangle_queue;

  destroyVertexBuffer();
  destroyTriangleBuffer();
}

static int sort_array_size = 0;
static SMvertex** sort_array = 0;

static void quicksort_for_border_markers(SMvertex** a, int i, int j)
{
  int in_i = i;
  int in_j = j;
  SMvertex* w;
  SMvertex* key = a[(i+j)/2];
  do
  {
    while ( a[i] < key ) i++;
    while ( a[j] > key ) j--;
    if (i<j)
    {
      w = a[i];
      a[i] = a[j];
      a[j] = w;
    }
  } while (++i<=--j);
  if (i == j+3)
  {
    i--;
    j++;
  }
  if (j>in_i) quicksort_for_border_markers(a, in_i, j);
  if (i<in_j) quicksort_for_border_markers(a, i, in_j);
}

static void finalize_vertex(SMvertex* vertex)
{
  int i,j,k;
  SMtriangle* triangle;
  // only vertices that have not yet been marked as boundary need checking
  if (vertex->boundary == false)
  {
    if (sort_array_size <= 2*vertex->list_size)
    {
      sort_array_size = 2*(vertex->list_size + 10);
      if (sort_array) free (sort_array);
      sort_array = (SMvertex**)malloc(sort_array_size*sizeof(SMvertex*));
    }
    k = 0;
    for (i = 0; i < vertex->list_size; i++)
    {
      triangle = vertex->list[i];
      j = 0;
      while (triangle->vertices[j] != vertex) j++;
      sort_array[k] = triangle->vertices[(j+1)%3];
      k++;
      sort_array[k] = triangle->vertices[(j+2)%3];
      k++;
    }
    // set the one past the last entry to zero to simplify the check below
    sort_array[k] = 0;
    // sort the array
    quicksort_for_border_markers(sort_array, 0, k-1);
    // check for boundary edges
    for (i = 0; i < k; i++)
    {
      if (sort_array[i] != sort_array[i+1])
      {
        vertex->boundary = true;
        sort_array[i]->boundary = true;
      }
      else
      {
        i++;
      }
    }
  }
  // check whether this makes any triangles ready for output
  for (i = 0; i < vertex->list_size; i++)
  {
    triangle = vertex->list[i];
    if (triangle->ready == 2)
    {
      triangle_queue->push_back(triangle);
    }
    else
    {
      triangle->ready++;
    }
  }
}

int SMreadWithBoundaryMarkers::read_input()
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
    // process the triangle's vertices
    for (i = 0; i < 3; i++)
    {
      // look for vertex in the hash
      hash_element = vertex_hash->find(smreader->t_idx[i]);
      // if vertex not in hash we have something other than a pre-order mesh  
      if (hash_element == vertex_hash->end())
      {
        // create a vertex and insert into hash
        vertex = allocVertex();
        // insert vertex into hash
        vertex_hash->insert(my_hash::value_type(smreader->t_idx[i], vertex));
        // increment the counter for the number vertices that are not in pre-order
        num_non_pre_order++;
      }
      else
      {
        // use vertex found in hash
        vertex = (*hash_element).second;
        // maybe we should remove it
        if (smreader->t_final[i]) vertex_hash->erase(hash_element);
      }
      // add vertex to triangle
      triangle->vertices[i] = vertex;
      // add triangle to vertex
      if (vertex->list_size == vertex->list_alloc)
      {
        vertex->list_alloc += 4;
        vertex->list = (SMtriangle**)realloc(vertex->list, sizeof(SMtriangle*)*vertex->list_alloc);
        if (!vertex->list)
        {
          fprintf(stderr, "FATAL ERROR: realloc for %d elements of vertex->list failed\n", vertex->list_alloc);
          return 0;
        }
      }
      vertex->list[vertex->list_size] = triangle;
      vertex->list_size++;
    }
    //  // does this triangle finalize any vertices
    for (i = 0; i < 3; i++)
    {
      // does this triangle finalize this vertex
      if (smreader->t_final[i])
      {
        // finalize the vertex
        finalize_vertex(triangle->vertices[i]);
      }
    }
  }
  else if (event == SM_VERTEX)
  {
    if (num_non_pre_order)
    {
      // look for vertex in the hash
      hash_element = vertex_hash->find(smreader->v_idx);
      if (hash_element == vertex_hash->end())
      {
        // create a vertex and insert into hash
        vertex = allocVertex();
        // insert vertex into hash
        vertex_hash->insert(my_hash::value_type(smreader->v_idx, vertex));
      }
      else
      {
        // decrement the counter for the number vertices that are not in pre-order
        num_non_pre_order--;
      }
      // copy vertex's coordinates
      VecCopy3fv(vertex->pos_f, smreader->v_pos_f);
    }
    else
    {
      // create a vertex and insert into hash
      vertex = allocVertex();
      // insert vertex into hash
      vertex_hash->insert(my_hash::value_type(smreader->v_idx, vertex));
      // copy vertex's coordinates
      VecCopy3fv(vertex->pos_f, smreader->v_pos_f);
    }
    // for post-order input this also finalizes the vertex
    if (smreader->post_order)
    {
      finalize_vertex(vertex);
    }
  }
  else if (event == SM_FINALIZED)
  {
    // look for finalized vertex in the vertex hash
    hash_element = vertex_hash->find(smreader->final_idx);
    // vertices must preceed their finalization
    if (hash_element == vertex_hash->end())
    {
      // fatal error
      fprintf(stderr, "FATAL ERROR: explicitely finalized vertex %d not in hash.\n", smreader->final_idx);
      return 0;
    }
    // get vertex found in hash
    vertex = (*hash_element).second;
    // erase its entry from the hash
    vertex_hash->erase(hash_element);
    // and finalize it
    finalize_vertex(vertex); 
    return 1;
  }
  else if (event == SM_EOF)
  {
    // all remaining vertices are implicitely finalized
    if (vertex_hash->size())
    {
      // get some hash element ... 
      hash_element = vertex_hash->begin();
      // get the vertex
      vertex = (*hash_element).second;
      // erase its entry from the hash
      vertex_hash->erase(hash_element);
      // and finalize it
      finalize_vertex(vertex); 
    }
    else
    {
      return 0;
    }
  }
  return 1;
}

int SMreadWithBoundaryMarkers::read_triangle()
{
  int i;
  SMvertex* vertex;
  SMtriangle* triangle;

  // read from input until tetrahedra are ready for ouput

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
	//    assert(vertex->index >= 0);
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
    // set additional vertex property array
    t_pos_f[i] = vertex->pos_f;
    // fill boundary marker array
    t_boundary[i] = vertex->boundary;
  }

  // deallocate this triangle

  deallocTriangle(triangle); // the memory can still be used until the next triangle is allocated.

  return 1;
}

SMevent SMreadWithBoundaryMarkers::read_element()
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
    // reference vertex properties in API interface
    v_idx = t_idx[new_vertices[next_new]];
    VecCopy3fv(v_pos_f,t_pos_f[new_vertices[next_new]]);
    v_boundary = t_boundary[new_vertices[next_new]];
    have_new--; next_new++;
    v_count++;
    return SM_VERTEX;
  }
  else
  {
    // reference triangle properties in API interface
    have_triangle = 0;
    f_count++;
    return SM_TRIANGLE;
  }
}

SMevent SMreadWithBoundaryMarkers::read_event()
{
  if ((have_new + have_triangle + have_finalized) == 0)
  {
    return read_element();
  }
  if (have_new)
  {
    // reference vertex properties in API interface
    v_idx = t_idx[new_vertices[next_new]];
    VecCopy3fv(v_pos_f,t_pos_f[new_vertices[next_new]]);
    v_boundary = t_boundary[new_vertices[next_new]];
    have_new--; next_new++;
    v_count++;
    return SM_VERTEX;
  }
  else if (have_triangle)
  {
    // reference triangle properties in API interface
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

SMreadWithBoundaryMarkers::SMreadWithBoundaryMarkers()
{
  // init of SMreadWithBoundaryMarkers
  v_boundary = false;
  t_pos_f[0] = t_pos_f[1] = t_pos_f[2] = 0;
  t_boundary[0] = t_boundary[1] = t_boundary[2] = false;

  smreader = 0;

  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;
}

SMreadWithBoundaryMarkers::~SMreadWithBoundaryMarkers()
{
  // clean-up for SMreadWithBoundaryMarkers
  bb_min_f = 0;
  bb_max_f = 0;
  if (smreader) delete smreader;
}
