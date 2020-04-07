/*
===============================================================================

  FILE:  smwriter_smd.cpp
  
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
#include "smwriter_smd.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "dynamicvector.h"
#include "dynamicqueue.h"
#include "littlecache.h"

#include "rangemodel.h"
#include "rangeencoder.h"

#include "floatcompressor.h"
#include "integercompressor_new.h"
#include "positionquantizer_new.h"

#include "vec3fv.h"
#include "vec3iv.h"

#include "hash_map.h"

#define PRINT_DEBUG_OUTPUT if (0) fprintf
#define PRINT_CONTROL_OUTPUT
#undef PRINT_CONTROL_OUTPUT

#define SMC_ADD 0
#define SMC_JOIN 1
#define SMC_FILL_END 2
#define SMC_SKIP 3
#define SMC_BORDER 4
#define SMC_START 3

// data structures used for streaming compression

struct SMedge;
struct SMtriangle;

typedef struct SMvertex
{
  SMvertex* buffer_next;    // used for efficient memory management (and also for the dynamic vector)
  // geometry
  int index;
  int use_total;
  float v[3];
  // incoming triangles
  int incoming_size;
  int incoming_alloc;
  SMtriangle** incoming;
  // used for compression
  int use_count;
  int list_size;
  int list_alloc;
  SMedge** list;
} SMvertex;

typedef struct SMtriangle
{
  SMtriangle* buffer_next;    // used for efficient memory management
  SMvertex* vertices[3];
} SMtriangle;

typedef struct SMedge
{
  SMedge* buffer_next;    // used for efficient memory management
  SMvertex* origin;
  SMvertex* target;
  float across[3];
} SMedge;

typedef hash_map<int, SMvertex*> my_hash;

static my_hash* vertex_hash;

static int next_waiting = 100;
static DynamicQueue* waiting_queue; // for triangles ready to be encoded
static DynamicQueue* traversal_queue; // for edges on the traversal front
static LittleCache* little_cache; // for subsequent traversal front misses?

static DynamicVector* dv;

static PositionQuantizerNew* pq;
static IntegerCompressorNew* ic[3];

static FloatCompressor* fc[3];

// rangecoder and probability tables

#define MAX_DEGREE_ONE 4
#define MAX_USE_COUNT 15
#define MAX_USE_COUNT_OP 10

static RangeEncoder* re_conn;
static RangeEncoder* re_conn_op;
static RangeEncoder* re_conn_rl;
static RangeEncoder* re_conn_index;
static RangeEncoder* re_conn_final;

static RangeEncoder* re_geom;

// is there more to encode
static RangeModel* rmDone;

// which operation
static RangeModel* rmWaitingOp;
static RangeModel*** rmTraversalOp;

// used for start operations
static RangeModel* rmWhichStart;

// used for fill/end operations
static RangeModel* rmRight;
static RangeModel* rmLeft;

// used for vertex finalization
static RangeModel*** rmFinalized;

static void initEncoder(FILE* file)
{
  if (file)
  {
    re_conn = new RangeEncoder(file);
    re_conn_rl = re_conn;
    re_conn_op = re_conn;
    re_conn_index = re_conn;
    re_conn_final = re_conn;
    re_geom = re_conn;
  }
  else
  {
    re_conn = new RangeEncoder(0,false);
    re_conn_op = new RangeEncoder(0,false);
    re_conn_rl = new RangeEncoder(0,false);
    re_conn_index = new RangeEncoder(0,false);
    re_conn_final = new RangeEncoder(0,false);
    re_geom = new RangeEncoder(0,false);
  }
}

static void finishEncoder(int nverts)
{
  if (re_conn != re_conn_op)
  {
    re_conn_op->done();
    re_conn_rl->done();
    re_conn_index->done();
    re_conn_final->done();
    re_conn->done();
    re_geom->done();

    fprintf(stderr,"\n");
    fprintf(stderr,"*** compression rate details ***\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"op:   \t%6.3f bpv\n", 8.0f/nverts*re_conn_op->getNumberBytes());
    fprintf(stderr,"rl:   \t%6.3f bpv\n", 8.0f/nverts*re_conn_rl->getNumberBytes());
    fprintf(stderr,"index:\t%6.3f bpv\n", 8.0f/nverts*re_conn_index->getNumberBytes());
    fprintf(stderr,"final:\t%6.3f bpv\n", 8.0f/nverts*re_conn_final->getNumberBytes());
    fprintf(stderr,"rest: \t%6.3f bpv\n", 8.0f/nverts*re_conn->getNumberBytes());
    fprintf(stderr,"\n");
    fprintf(stderr,"conn: \t%6.3f bpv\n", 8.0f/nverts*(re_conn_op->getNumberBytes()+re_conn_rl->getNumberBytes()+re_conn_index->getNumberBytes()+re_conn_final->getNumberBytes()+re_conn->getNumberBytes()));
    fprintf(stderr,"geom: \t%6.3f bpv\n", 8.0f/nverts*re_geom->getNumberBytes());
    fprintf(stderr,"\n");
    fprintf(stderr,"total:\t%6.3f bpv\n", 8.0f/nverts*(re_geom->getNumberBytes()+re_conn_op->getNumberBytes()+re_conn_rl->getNumberBytes()+re_conn_index->getNumberBytes()+re_conn_final->getNumberBytes()+re_conn->getNumberBytes()));

    delete re_conn_op;
    delete re_conn_rl;
    delete re_conn_index;
    delete re_conn_final;
    delete re_conn;
    delete re_geom;
  }
  else
  {
    re_conn->done();
#ifdef PRINT_CONTROL_OUTPUT
    fprintf(stderr,"total: bytes %d bpv %f\n", re_conn->getNumberBytes(),(float)re_conn->getNumberBits()/nverts);
#endif
    delete re_conn;
  }
}

static void initModels(int compress)
{
  int i,j;

  // range table for whether there is more
  rmDone = new RangeModel(2,0,compress);

  // range table for which operation it is
  rmWaitingOp = new RangeModel(4,0,compress);
  rmTraversalOp = (RangeModel***)malloc(sizeof(RangeModel**)*MAX_USE_COUNT_OP);
  for (i = 1; i < MAX_USE_COUNT_OP; i++)
  {
    rmTraversalOp[i] = (RangeModel**)malloc(sizeof(RangeModel*)*MAX_USE_COUNT_OP);
    for (j = 1; j < MAX_USE_COUNT_OP; j++)
    {
      rmTraversalOp[i][j] = new RangeModel(5,0,compress);
    }
  }

  // range tables for which start operation it is
  rmWhichStart = new RangeModel(4,0,compress);

  // range tables for fille (F) / end (E) operations
  rmRight = new RangeModel(2,0,compress);
  rmLeft = new RangeModel(2,0,compress);

  rmFinalized = (RangeModel***)malloc(sizeof(RangeModel**)*MAX_DEGREE_ONE);
  for (i = 0; i < MAX_DEGREE_ONE; i++)
  {
    rmFinalized[i] = (RangeModel**)malloc(sizeof(RangeModel*)*MAX_USE_COUNT);
    for (j = 1; j < MAX_USE_COUNT; j++)
    {
      rmFinalized[i][j] = new RangeModel(2,0,compress);
    }
  }
}

static void finishModels()
{
  int i,j;

  delete rmDone;

  // range table for which operation it is
  delete rmWaitingOp;
  for (i = 1; i < MAX_USE_COUNT_OP; i++)
  {
    for (j = 1; j < MAX_USE_COUNT_OP; j++)
    {
      delete rmTraversalOp[i][j];
    }
    free(rmTraversalOp[i]);
  }
  free(rmTraversalOp);

  delete rmWhichStart;

  delete rmRight;
  delete rmLeft;

  for (i = 0; i < MAX_DEGREE_ONE; i++)
  {
    for (j = 1; j < MAX_USE_COUNT; j++)
    {
      delete rmFinalized[i][j];
    }
    free(rmFinalized[i]);
  }
  free(rmFinalized);
}

// statistics

#ifdef PRINT_CONTROL_OUTPUT
static int op_start = 0;
static int op_add = 0;
static int op_join = 0;
static int op_fill = 0;
static int op_end = 0;
static int op_skip = 0;
static int op_border = 0;

static int used_index = 0;

static int prediction_none = 0;
static int prediction_last = 0;
static int prediction_across = 0;

static int right_confirm = 0;
static int right_correct = 0;
static int left_confirm = 0;
static int left_correct = 0;

static int max_in_width = 0;
static int max_in_span = 0;
static int max_out_width = 0;
static int max_out_span = 0;
static int v_out_count = 0;

static int vertex_buffer_maxsize = 0;
static int edge_buffer_maxsize = 0;
static int triangle_buffer_maxsize = 0;
#endif

// efficient memory allocation

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
    vertex_buffer_next[i].list = 0;
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
      vertex_buffer_next[i].list = 0;
    }
    vertex_buffer_next[vertex_buffer_alloc-1].buffer_next = 0;
    vertex_buffer_alloc = 2*vertex_buffer_alloc;
  }
  // get pointer to next available vertex
  SMvertex* vertex = vertex_buffer_next;
  vertex_buffer_next = vertex->buffer_next;

  vertex->use_total = 0;
  if (vertex->incoming == 0)
  {
    vertex->incoming = (SMtriangle**)malloc(sizeof(SMtriangle*)*10);
    vertex->incoming_alloc = 10;
  }
  vertex->incoming_size = 0;
  vertex->use_count = 0;
  if (vertex->list == 0)
  {
    vertex->list = (SMedge**)malloc(sizeof(SMedge*)*6);
    vertex->list_alloc = 6;
  }
  vertex->list_size = 0;

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

static int edge_buffer_size = 0;
static int edge_buffer_alloc = 1024;
static SMedge* edge_buffer_next = 0;

static int initEdgeBuffer(int size)
{
  edge_buffer_next = (SMedge*)malloc(sizeof(SMedge)*size);

  if (edge_buffer_next == 0)
  {
    fprintf(stderr,"malloc for edge buffer failed\n");
    return 0;
  }
  for (int i = 0; i < size; i++)
  {
    edge_buffer_next[i].buffer_next = &(edge_buffer_next[i+1]);
  }
  edge_buffer_next[size-1].buffer_next = 0;
  edge_buffer_alloc = size;
  edge_buffer_size = 0;
  return 1;
}

static SMedge* allocEdge(float* v)
{
  if (edge_buffer_next == 0)
  {
    edge_buffer_next = (SMedge*)malloc(sizeof(SMedge)*edge_buffer_alloc);
    if (edge_buffer_next == 0)
    {
      fprintf(stderr,"malloc for edge buffer failed\n");
      return 0;
    }
    for (int i = 0; i < edge_buffer_alloc; i++)
    {
      edge_buffer_next[i].buffer_next = &(edge_buffer_next[i+1]);
    }
    edge_buffer_next[edge_buffer_alloc-1].buffer_next = 0;
    edge_buffer_alloc = 2*edge_buffer_alloc;
  }
  // get index of next available vertex
  SMedge* edge = edge_buffer_next;
  edge_buffer_next = edge->buffer_next;
 
  edge->origin = 0;
  edge->target = 0;
  VecCopy3fv(edge->across,v);

  edge_buffer_size++;
  
#ifdef PRINT_CONTROL_OUTPUT
  if (edge_buffer_size > edge_buffer_maxsize) edge_buffer_maxsize = edge_buffer_size;
#endif

  return edge;
}

static void deallocEdge(SMedge* edge)
{
  edge->buffer_next = edge_buffer_next;
  edge_buffer_next = edge;
  edge_buffer_size--;
}

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

// helper functions

static void compressVertexPosition(float* n)
{
#ifdef PRINT_CONTROL_OUTPUT
  prediction_none++;
#endif

  if (pq)
  {
    pq->EnQuantize(n, (int*)n);
    for (int i = 0; i < 3; i++)
    {
      ic[i]->CompressNone(((int*)n)[i]);
    }
  }
  else
  {
    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->CompressNone(n[i]);
    }
  }
}

static void compressVertexPosition(float* l, float* n)
{
#ifdef PRINT_CONTROL_OUTPUT
  prediction_last++;
#endif

  if (pq)
  {
    pq->EnQuantize(n, (int*)n);
    for (int i = 0; i < 3; i++)
    {
      ic[i]->CompressLast(((int*)l)[i],((int*)n)[i]);
    }
  }
  else
  {
    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->CompressLast(l[i],n[i]);
    }
  }
}

static void compressVertexPosition(const float* a, const float* b, const float* c, float* n)
{
#ifdef PRINT_CONTROL_OUTPUT
  prediction_across++;
#endif

  if (pq)
  {
    int pred[3];
    VecAdd3iv(pred, (const int*)a, (const int*)c);
    VecSelfSubtract3iv(pred, (const int*)b);
    pq->Clamp(pred);
    pq->EnQuantize(n, (int*)n);
    for (int i = 0; i < 3; i++)
    {
      ic[i]->CompressAcross(pred[i], ((int*)n)[i]);
    }
  }
  else
  {
    float pred[3];
    VecParallelogram3fv(pred, a, b, c);

    /*
    This is CRAZY Micosoft crap.... the floating-point calculation below can have different results in Debug and Release mode!!!!
    VecAdd3fv(pred, a, c);
    VecSelfSubtract3fv(pred, b);
    */

    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->CompressAcross(pred[i],n[i]);
    }
  }
}

static void addEdgeToVertex(SMedge* edge, SMvertex* vertex)
{
  if (vertex->list_size == vertex->list_alloc)
  {
    vertex->list_alloc += 2;
    vertex->list = (SMedge**)realloc(vertex->list, sizeof(SMedge*)*vertex->list_alloc);
  }
  vertex->list[vertex->list_size] = edge;
  vertex->list_size++;
}

static void removeEdgeFromVertex(SMedge* edge, SMvertex* vertex)
{
  for (int i = 0; i < vertex->list_size; i++)
  {
    if (vertex->list[i] == edge)
    {
      vertex->list_size--;
      vertex->list[i] = vertex->list[vertex->list_size];
      return;
    }
  }
  fprintf(stderr,"ERROR: did not find edge in vertex list\n");
  exit(0);
}

static void removeEdgeFromVertices(SMedge* edge)
{
  removeEdgeFromVertex(edge, edge->origin);
  removeEdgeFromVertex(edge, edge->target);
}

static void addTriangleToVertex(SMtriangle* triangle, SMvertex* vertex)
{
  if (vertex->incoming_size == vertex->incoming_alloc)
  {
    vertex->incoming_alloc += 2;
    vertex->incoming = (SMtriangle**)realloc(vertex->incoming, sizeof(SMtriangle*)*vertex->incoming_alloc);
  }
  vertex->incoming[vertex->incoming_size] = triangle;
  vertex->incoming_size++;
  vertex->use_total--;
}

static void removeTriangleFromVertex(SMtriangle* triangle, SMvertex* vertex)
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
  fprintf(stderr,"ERROR: did not find triangle in vertex list\n");
  exit(0);
}

void SMwriter_smd::write_vertex(const float* v_pos_f)
{
  SMvertex* vertex = allocVertex();
  VecCopy3fv(vertex->v, v_pos_f);
  vertex_hash->insert(my_hash::value_type(v_count, vertex));
  v_count++;
#ifdef PRINT_CONTROL_OUTPUT
  if (vertex_hash->size() > max_in_width) max_in_width = vertex_hash->size();
#endif
}

bool SMwriter_smd::compress_triangle_waiting()
{
  int i, j;
  int rot = 0;
  int candidate, candidate_count;
  int degree_one, use_count;

  SMedge* edge;
  SMedge* edges[3];
  SMvertex* vertex;
  SMvertex* vertices[3];

  SMtriangle* triangle = (SMtriangle*)waiting_queue->getAndRemoveFirstElement();

  // are vertices already visited
  int num_v_visited = 0;
  int last_v_visited;
  int last_v_unvisited;

  for (i = 0; i < 3; i++)
  {
    if (triangle->vertices[i]->use_count)
    {
      vertices[i] = triangle->vertices[i];
      last_v_visited = i;
      num_v_visited++;
    }
    else
    {
      vertices[i] = 0;
      last_v_unvisited = i;
    }
  }

  // are edges already visited
  SMedge* e_visited[3];
  int num_e_visited = 0;
  int last_e_visited = -1;
  int last_e_unvisited;

  for (i = 0; i < 3; i++)
  {
    e_visited[i] = 0;
    if (vertices[i] && vertices[(i+1)%3]) // check if an inverse edge already exists
    {
      for (j = 0; j < vertices[i]->list_size; j++)
      {
        edge = vertices[i]->list[j];
        if (edge->origin == vertices[(i+1)%3])
        {
          e_visited[i] = edge;
          last_e_visited = i;
          num_e_visited++;
          break;
        }
      }
      if (last_e_visited != i)
      {
        last_e_unvisited = i;
      }
    }
  }

  // compress triangle depending on its adjacency

  if (num_e_visited == 0) // start
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_start++;
#endif

    if (dv->size() == 0)
    {
      re_conn_op->encode(rmDone, 0); // not done
    }
    else
    {
      re_conn_op->encode(rmWaitingOp, SMC_START);
    }

    PRINT_DEBUG_OUTPUT(stderr, "w%d ", SMC_START);

    // which START_i operation is it (0 <= i <= 3) 

    re_conn_op->encode(rmWhichStart, num_v_visited);

    // rotate triangle if necessary

    if (num_v_visited == 1)
    {
      rot = last_v_visited;
    }
    else if (num_v_visited == 2)
    {
      rot = (last_v_unvisited + 1) % 3;
    }

    // encode the (rotated) start configuration

    vertices[0] = triangle->vertices[rot];
    vertices[1] = triangle->vertices[(rot+1)%3];
    vertices[2] = triangle->vertices[(rot+2)%3];

    edges[0] = 0;
    edges[1] = 0;
    edges[2] = 0;

    // encode vertex v0

    if (vertices[0]->use_count)
    {
      // encode its index among all active vertices
      i = dv->getRelativeIndex(vertices[0]);
      re_conn_index->encode(dv->size(),i);
#ifdef PRINT_CONTROL_OUTPUT
      used_index++;
#endif
      PRINT_DEBUG_OUTPUT(stderr, "SNM0 ");
    }
    else
    {
      // encode its position
      compressVertexPosition(vertices[0]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[0]);
#ifdef PRINT_CONTROL_OUTPUT
      if (dv->size() > max_out_width) max_out_width = dv->size();
      vertices[0]->index = v_out_count;
      v_out_count++;
#endif
    }

    // encode vertex v1

    if (vertices[1]->use_count)
    {
      // encode its index among all active vertices
      i = dv->getRelativeIndex(vertices[1]);
      re_conn_index->encode(dv->size(),i);        
#ifdef PRINT_CONTROL_OUTPUT
      used_index++;
#endif
      PRINT_DEBUG_OUTPUT(stderr, "SNM1 ");
    }
    else
    {
      // encode its position
      compressVertexPosition(vertices[0]->v, vertices[1]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[1]);
#ifdef PRINT_CONTROL_OUTPUT
      if (dv->size() > max_out_width) max_out_width = dv->size();
      vertices[1]->index = v_out_count;
      v_out_count++;
#endif
    }

    // encode vertex v2

    if (vertices[2]->use_count)
    {
      // encode its index among all active vertices
      i = dv->getRelativeIndex(vertices[2]);
      re_conn_index->encode(dv->size(),i);        
#ifdef PRINT_CONTROL_OUTPUT
      used_index++;
#endif
      PRINT_DEBUG_OUTPUT(stderr, "SNM2 ");
    }
    else
    {
      // encode its position
      compressVertexPosition(vertices[0]->v, vertices[2]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[2]);
#ifdef PRINT_CONTROL_OUTPUT
      if (dv->size() > max_out_width) max_out_width = dv->size();
      vertices[2]->index = v_out_count;
      v_out_count++;
#endif
    }
  }
  else if (num_e_visited == 1) // add or join
  {
    if (num_v_visited == 2) // add
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_add++;
#endif
      re_conn_op->encode(rmWaitingOp, SMC_ADD);
      PRINT_DEBUG_OUTPUT(stderr, "w%d ", SMC_ADD);
    }
    else
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_join++;
#endif
      re_conn_op->encode(rmWaitingOp, SMC_JOIN);
      PRINT_DEBUG_OUTPUT(stderr, "w%d ", SMC_JOIN);
    }

    // rotate triangle if necessary

    rot = last_e_visited;

    // encode the (rotated) add (or join) configuration

    vertices[0] = triangle->vertices[rot];
    vertices[1] = triangle->vertices[(rot+1)%3];
    vertices[2] = triangle->vertices[(rot+2)%3];

    edges[0] = e_visited[rot];
    edges[1] = 0;
    edges[2] = 0;

    // encode the index of v0
    i = dv->getRelativeIndex(vertices[0]);
    re_conn_index->encode(dv->size(),i);        
#ifdef PRINT_CONTROL_OUTPUT
    used_index++;
#endif
    PRINT_DEBUG_OUTPUT(stderr, "m%d ", i);

    // encode which of the edges incident to v0 is e0
    candidate_count = 0;

    // how many potential candidates are around that vertex?
    for (i = 0; i < vertices[0]->list_size; i++)
    {
      if (vertices[0]->list[i]->target == vertices[0])
      {
        if (vertices[0]->list[i] == edges[0]) candidate = candidate_count;
        candidate_count++;
      }
    }

    // only if there is more than one we need to encode it
    if (candidate_count > 1)
    {
      re_conn->encode(candidate_count,candidate);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }

    // encode v2
    if (vertices[2]->use_count) // join
    {
      // encode its index among all active vertices
      i = dv->getRelativeIndex(vertices[2]);
      re_conn_index->encode(dv->size(),i);        
#ifdef PRINT_CONTROL_OUTPUT
      used_index++;
#endif
      PRINT_DEBUG_OUTPUT(stderr, "J %d ",i);
    }
    else // add
    {
      // encode its position
      compressVertexPosition(vertices[0]->v, edges[0]->across, vertices[1]->v, vertices[2]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[2]);
#ifdef PRINT_CONTROL_OUTPUT
      if (dv->size() > max_out_width) max_out_width = dv->size();
      vertices[2]->index = v_out_count;
      v_out_count++;
#endif
    }
  }
  else // fill or end
  {
    re_conn_op->encode(rmWaitingOp, SMC_FILL_END);

    PRINT_DEBUG_OUTPUT(stderr, "w%d ", SMC_FILL_END);

    // rotate triangle if necessary

    if (num_e_visited == 2)
    {
      rot = (last_e_unvisited + 2) % 3;
    }

    // encode the (rotated) fill (or end) configuration

    vertices[0] = triangle->vertices[rot];
    vertices[1] = triangle->vertices[(rot+1)%3];
    vertices[2] = triangle->vertices[(rot+2)%3];

    edges[0] = e_visited[rot];
    edges[1] = e_visited[(rot+1)%3];
    edges[2] = e_visited[(rot+2)%3];

    // encode the index of v0
    i = dv->getRelativeIndex(vertices[0]);
    re_conn_index->encode(dv->size(),i);
#ifdef PRINT_CONTROL_OUTPUT
    used_index++;
#endif
    PRINT_DEBUG_OUTPUT(stderr, "m%d ", i);

    // encode edge e0 (e.g. the edge from v0 to v1)
    candidate_count = 0;

    // how many potential candidates are around that vertex?
    for (i = 0; i < vertices[0]->list_size; i++)
    {
      if (vertices[0]->list[i]->target == vertices[0])
      {
        if (vertices[0]->list[i] == edges[0]) candidate = candidate_count;
        candidate_count++;
      }
    }

    // only if there is more than one we need to encode it
    if (candidate_count > 1)
    {
      re_conn->encode(candidate_count,candidate);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }

    // encode edge e2 (e.g. the edge from v2 to v0)
    candidate_count = 0;

    // how many potential candidates are around that vertex?
    for (i = 0; i < vertices[0]->list_size; i++)
    {
      if (vertices[0]->list[i]->origin == vertices[0])
      {
        if (vertices[0]->list[i] == edges[2]) candidate = candidate_count;
        candidate_count++;
      }
    }

    // only if there is more than one we need to encode it
    if (candidate_count > 1)
    {
      re_conn->encode(candidate_count,candidate);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }

#ifdef PRINT_CONTROL_OUTPUT
    if (num_e_visited == 2)
    {
      op_fill++;
    }
    else
    {
      op_end++;
    }
#endif
  }

  // increment vertex use_counts, create edges, and update edge degrees

  for (i = 0; i < 3; i++)
  {
    vertices[i]->use_count++;

    if (edges[i] == 0)
    {
      edges[i] = allocEdge(vertices[(i+2)%3]->v); // the vertex coordinates are used for parallelogram predictions
      edges[i]->origin = vertices[i];
      edges[i]->target = vertices[(i+1)%3];
      addEdgeToVertex(edges[i],vertices[i]);
      addEdgeToVertex(edges[i],vertices[(i+1)%3]);
      traversal_queue->addElement(edges[i]);
    }
    else
    {
      if (((int*)edges[i])[0] >= 0)
      {
        traversal_queue->removeElement(edges[i]);
      }
      removeEdgeFromVertices(edges[i]);
      deallocEdge(edges[i]);
    }
  }

  // write finalization info, finalize out data structures, delete triangle pointer

  for (i = 0; i < 3; i++)
  {
    vertex = vertices[i];

    degree_one = vertex->list_size;
    if (degree_one >= MAX_DEGREE_ONE)
    {
      degree_one = MAX_DEGREE_ONE-1;
    }
    use_count = vertex->use_count;
    if (use_count >= MAX_USE_COUNT)
    {
      use_count = MAX_USE_COUNT-1;
    }

    if (vertex->use_count == vertex->use_total) // does this triangle finalize this vertex?
    {
      re_conn_final->encode(rmFinalized[degree_one][use_count], 1);
      dv->removeElement(vertex);
      for (j = 0; j < vertex->list_size; j++)
      {
        edge = vertex->list[j];
        if (vertex == edge->origin)
        {
          removeEdgeFromVertex(edge, edge->target);
        }
        else
        {
          removeEdgeFromVertex(edge, edge->origin);
        }
        if (((int*)edge)[0] >= 0) traversal_queue->removeElement(edge);
        deallocEdge(edge);
      }
#ifdef PRINT_CONTROL_OUTPUT
      if ((v_out_count-vertex->index+1) > max_out_span) max_out_span = (v_out_count-vertex->index+1);
#endif
      deallocVertex(vertex);
    }
    else
    {
      re_conn_final->encode(rmFinalized[degree_one][use_count], 0);
      removeTriangleFromVertex(triangle, vertex);
    }
  }
  deallocTriangle(triangle);
  return true;
}

bool SMwriter_smd::compress_triangle_traversal()
{
  int i,j;
  int candidate,candidate_count;

  SMedge* edges[3];
  SMvertex* vertices[3];

  SMtriangle* triangle;
  SMvertex* vertex;
  SMedge* edge = (SMedge*)traversal_queue->getAndRemoveFirstElement();

  // check if there is a triangle for us to compress
  vertices[0] = edge->target;
  vertices[1] = edge->origin;
  vertices[2] = 0;
  for (i = 0; i < vertices[0]->incoming_size; i++)
  {
    triangle = vertices[0]->incoming[i];
    if (triangle->vertices[0] == vertices[0] && triangle->vertices[1] == vertices[1])
    {
      vertices[2] = triangle->vertices[2];
      break;
    }
    if (triangle->vertices[1] == vertices[0] && triangle->vertices[2] == vertices[1])
    {
      vertices[2] = triangle->vertices[0];
      break;
    }
    if (triangle->vertices[2] == vertices[0] && triangle->vertices[0] == vertices[1])
    {
      vertices[2] = triangle->vertices[1];
      break;
    }
  }
  if (vertices[2] == 0) // there is no triangle adjacent to this edge. why?
  {
    if ((vertices[0]->use_total > 0) && (vertices[1]->use_total > 0)) // because it is a border edge
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_border++;
#endif
      i = (vertices[0]->use_count < MAX_USE_COUNT_OP)?vertices[0]->use_count:MAX_USE_COUNT_OP-1;
      j = (vertices[1]->use_count < MAX_USE_COUNT_OP)?vertices[1]->use_count:MAX_USE_COUNT_OP-1;
      re_conn_op->encode(rmTraversalOp[i][j], SMC_BORDER); 
      removeEdgeFromVertices(edge);
      deallocEdge(edge);
      PRINT_DEBUG_OUTPUT(stderr, "t%d:%d:%d ",i,j,SMC_BORDER);
      if (traversal_queue->elements())
      {
        return compress_triangle_traversal();
      }
      return false;
    }
    else // because it may come later 
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_skip++;
#endif
      i = (vertices[0]->use_count < MAX_USE_COUNT_OP)?vertices[0]->use_count:MAX_USE_COUNT_OP-1;
      j = (vertices[1]->use_count < MAX_USE_COUNT_OP)?vertices[1]->use_count:MAX_USE_COUNT_OP-1;
      re_conn_op->encode(rmTraversalOp[i][j], SMC_SKIP); 
      ((int*)edge)[0] = -1; // mark element as no longer being in the queue
      PRINT_DEBUG_OUTPUT(stderr, "t%d:%d:%d ",i,j,SMC_SKIP);
      return false;
    }
  }
  
  // is this triangle adjacent to any other edge?

  edges[0] = edge;
  edges[1] = 0;
  edges[2] = 0;

  // check by looping around v2

  for (i = 0; i < vertices[2]->list_size; i++)
  {
    edge = vertices[2]->list[i];
    if (vertices[0] == edge->origin)
    {
      edges[2] = edge;
    }
    else if (vertices[1] == edge->target)
    {
      edges[1] = edge;
    }
  }


  // compress the triangle based on its operation: add/join, or fill/end

  if (edges[1] || edges[2]) // fill or end
  {
    i = (vertices[0]->use_count < MAX_USE_COUNT_OP)?vertices[0]->use_count:MAX_USE_COUNT_OP-1;
    j = (vertices[1]->use_count < MAX_USE_COUNT_OP)?vertices[1]->use_count:MAX_USE_COUNT_OP-1;
    re_conn_op->encode(rmTraversalOp[i][j], SMC_FILL_END); 
    PRINT_DEBUG_OUTPUT(stderr, "t%d:%d:%d ",i,j,SMC_FILL_END);

    if (edges[1] && edges[2])
    {
      if (vertices[1]->use_count >= vertices[0]->use_count)
      {
        i = 1;
        PRINT_DEBUG_OUTPUT(stderr, "r0 ");
        re_conn_rl->encode(rmRight, 0); // confirm
#ifdef PRINT_CONTROL_OUTPUT
        right_confirm++;
#endif
      }
      else
      {
        i = 0;
        PRINT_DEBUG_OUTPUT(stderr, "l0 ");
        re_conn_rl->encode(rmLeft, 0); // confirm
#ifdef PRINT_CONTROL_OUTPUT
        left_confirm++;
#endif
      }
#ifdef PRINT_CONTROL_OUTPUT
      op_end++;
#endif
    }
    else
    {
      if (vertices[1]->use_count >= vertices[0]->use_count)
      {
        if (edges[1])
        {
          i = 1;
          PRINT_DEBUG_OUTPUT(stderr, "r0 ");
          re_conn_rl->encode(rmRight, 0); // confirm
#ifdef PRINT_CONTROL_OUTPUT
          right_confirm++;
#endif
        }
        else
        {
          i = 0;
          PRINT_DEBUG_OUTPUT(stderr, "r1 ");
          re_conn_rl->encode(rmRight, 1); // correct
#ifdef PRINT_CONTROL_OUTPUT
          right_correct++;
#endif
        }
      }
      else
      {
        if (edges[2])
        {
          i = 0;
          PRINT_DEBUG_OUTPUT(stderr, "l0 ");
          re_conn_rl->encode(rmLeft, 0); // confirm
#ifdef PRINT_CONTROL_OUTPUT
          left_confirm++;
#endif
        }
        else
        {
          i = 1;
          PRINT_DEBUG_OUTPUT(stderr, "l1 ");
          re_conn_rl->encode(rmLeft, 1); // correct
#ifdef PRINT_CONTROL_OUTPUT
          left_correct++;
#endif
        }
      }
#ifdef PRINT_CONTROL_OUTPUT
      op_fill++;
#endif
    }

    if (i == 1)
    {
      candidate_count = 0;

      // how many edges are around v1 that could potentially be e1?
      for (i = 0; i < vertices[1]->list_size; i++)
      {
        if (vertices[1]->list[i]->target == vertices[1])
        {
          if (vertices[1]->list[i] == edges[1]) candidate = candidate_count;
          candidate_count++;
        }
      }
      // only if there is more than one we need to encode it
      if (candidate_count > 1)
      {
        re_conn->encode(candidate_count,candidate);
        PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
      }
    }
    else
    {
      candidate_count = 0;
      
      // how many edges are around v0 that could potentially be e2?
      for (i = 0; i < vertices[0]->list_size; i++)
      {
        if (vertices[0]->list[i]->origin == vertices[0])
        {
          if (vertices[0]->list[i] == edges[2]) candidate = candidate_count;
          candidate_count++;
        }
      }

      // only if there is more than one we need to encode it
      if (candidate_count > 1)
      {
        re_conn->encode(candidate_count,candidate);
        PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
      }
    }
  }
  else
  {
    // encode v2
    if (vertices[2]->use_count)
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_join++;
#endif
      i = (vertices[0]->use_count < MAX_USE_COUNT_OP)?vertices[0]->use_count:MAX_USE_COUNT_OP-1;
      j = (vertices[1]->use_count < MAX_USE_COUNT_OP)?vertices[1]->use_count:MAX_USE_COUNT_OP-1;
      re_conn_op->encode(rmTraversalOp[i][j], SMC_JOIN); 
      PRINT_DEBUG_OUTPUT(stderr, "t%d:%d:%d ",i,j,SMC_JOIN);

      // get v2's index among all active vertices
      i = dv->getRelativeIndex(vertices[2]);
      // encode v2's index
      re_conn_index->encode(dv->size(),i);        
    }
    else
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_add++;
#endif
      i = (vertices[0]->use_count < MAX_USE_COUNT_OP)?vertices[0]->use_count:MAX_USE_COUNT_OP-1;
      j = (vertices[1]->use_count < MAX_USE_COUNT_OP)?vertices[1]->use_count:MAX_USE_COUNT_OP-1;
      re_conn_op->encode(rmTraversalOp[i][j], SMC_ADD); 
      PRINT_DEBUG_OUTPUT(stderr, "t%d:%d:%d ",i,j,SMC_ADD);

      // encode v2's position
      compressVertexPosition(vertices[0]->v, edges[0]->across, vertices[1]->v, vertices[2]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[2]);
#ifdef PRINT_CONTROL_OUTPUT
      if (dv->size() > max_out_width) max_out_width = dv->size();
      vertices[2]->index = v_out_count;
      v_out_count++;
#endif
    }
  }

  // either create edges and add to queue or delete edges and remove from queue

  removeEdgeFromVertices(edges[0]);
  deallocEdge(edges[0]);

  // note ... edges[0] is treated already as edge

  for (i = 1; i < 3; i++)
  {
    if (edges[i] == 0)
    { 
      edges[i] = allocEdge(vertices[(i+2)%3]->v);  // the vertex coordinates are used for parallelogram predictions
      edges[i]->origin = vertices[i];
      edges[i]->target = vertices[(i+1)%3];
      addEdgeToVertex(edges[i],vertices[i]);
      addEdgeToVertex(edges[i],vertices[(i+1)%3]);
      traversal_queue->addElement(edges[i]);
    }
    else
    {
      if (((int*)edges[i])[0] >= 0)
      {
        traversal_queue->removeElement(edges[i]);
        ((int*)edges[i])[0] = -1;
      }
      removeEdgeFromVertices(edges[i]);
      deallocEdge(edges[i]);
    }
  }

  // increment vertex use_counts, write finalization info and finalize data structures

  int degree_one;
  int use_count;

  for (i = 0; i < 3; i++)
  {
    vertex = vertices[i];
    vertex->use_count++;

    degree_one = vertex->list_size;
    if (degree_one >= MAX_DEGREE_ONE)
    {
      degree_one = MAX_DEGREE_ONE-1;
    }
    use_count = vertex->use_count;
    if (use_count >= MAX_USE_COUNT)
    {
      use_count = MAX_USE_COUNT-1;
    }
    
    if (vertex->use_count == vertex->use_total) // does this triangle finalize this vertex?
    {
      re_conn_final->encode(rmFinalized[degree_one][use_count], 1);
      dv->removeElement(vertex);
      for (j = 0; j < vertex->list_size; j++)
      {
        edge = vertex->list[j];
        if (vertex == edge->origin)
        {
          removeEdgeFromVertex(edge, edge->target);
        }
        else
        {
          removeEdgeFromVertex(edge, edge->origin);
        }
        if (((int*)edge)[0] >= 0) traversal_queue->removeElement(edge);
        deallocEdge(edge);
      }
#ifdef PRINT_CONTROL_OUTPUT
      if ((v_out_count-vertex->index+1) > max_out_span) max_out_span = (v_out_count-vertex->index+1);
#endif
      deallocVertex(vertex);
    }
    else
    {
      re_conn_final->encode(rmFinalized[degree_one][use_count], 0);
      removeTriangleFromVertex(triangle, vertex);
    }
  }
  waiting_queue->removeElement(triangle);
  deallocTriangle(triangle);
  return true;
}

bool SMwriter_smd::compress_triangle()
{
  if (traversal_queue->elements())
  {
    if (next_waiting > 0)
    {
      if (compress_triangle_traversal())
      {
        next_waiting--;
        return true;
      }
    }
  }

  if (waiting_queue->elements())
  {
    if (compress_triangle_waiting())
    {
      next_waiting = (100 < traversal_queue->elements() ? traversal_queue->elements() : 100);
      return true;
    }
  }

  return false;
}

void SMwriter_smd::write_triangle(const int* t_idx)
{
  bool t_final[] = {false, false, false};
  write_triangle(t_idx, t_final);
}

void SMwriter_smd::write_triangle(const int* t_idx, const bool* t_final)
{
  if (f_count == 0) write_header();

  int i;
  my_hash::iterator hash_elements[3];
  SMtriangle* triangle = allocTriangle();

  // get vertices from hash
  for (i = 0; i < 3; i++)
  {
    hash_elements[i] = vertex_hash->find(t_idx[i]);
    if (hash_elements[i] == vertex_hash->end())
    {
      triangle->vertices[i] = allocVertex();
      vertex_hash->insert(my_hash::value_type(t_idx[i], triangle->vertices[i]));
    }
    else
    {
      triangle->vertices[i] = (*hash_elements[i]).second;
      addTriangleToVertex(triangle,triangle->vertices[i]);
    }
  }

  // check for finalization
  for (i = 0; i < 3; i++)
  {
    if (t_final[i])
    {
      triangle->vertices[i]->use_total *= -1;
      vertex_hash->erase(hash_elements[i]);
#ifdef PRINT_CONTROL_OUTPUT
      if ((v_count - t_idx[i] + 1) > max_in_span) max_in_span = (v_count - t_idx[i] + 1);
#endif
    }
  }

  // add triangle to buffer
  waiting_queue->addElement(triangle);

  // compress triangles while the buffer is full
  if (max_delay > 0)
  {
    while (waiting_queue->size() >= max_delay)
    {
      compress_triangle();
    }
  }
  else
  {
    while (waiting_queue->size() > vertex_hash->size()*width_delay)
    {
      compress_triangle();
    }
  }

  f_count++;
}

void SMwriter_smd::write_finalized(int final_idx)
{
  my_hash::iterator hash_element = vertex_hash->find(final_idx);
  if (hash_element == vertex_hash->end())
  {
    fprintf(stderr,"FATAL ERROR: finalized vertex %d not in hash\n",final_idx);
    exit(0);
  }
  (*hash_element).second->use_total *= -1;
  vertex_hash->erase(hash_element);
}

#define SM_VERSION 2 // this is SMD

bool SMwriter_smd::open(FILE* file, int nbits, int delay)
{
#ifdef _WIN32
  if (file == stdout)
  {
    if(_setmode( _fileno( stdout ), _O_BINARY ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdout to binary (untranslated) mode\n");
    }
  }
#endif

  if (nbits < 0 || nbits > 24)
  {
    fprintf(stderr, "ERROR: nbits is %d. it should be in the range [0,24]\n");
    exit(0);
  }

  // write version
  if (file)
  {
    fputc(SM_VERSION, file);
  }

  initVertexBuffer(1024);
  initEdgeBuffer(1024);
  initTriangleBuffer(1024);

  waiting_queue = new DynamicQueue();
  traversal_queue = new DynamicQueue();
  little_cache = new LittleCache();

  dv = new DynamicVector();
  vertex_hash = new my_hash;

  initEncoder(file);
  initModels(1);

  // write precision
  re_conn->encode(25,nbits);

  // store precision (is needed in write_header)
  this->nbits = nbits;

  if (delay < 0)
  {
    max_delay = 0;
    width_delay = -delay;
    fprintf(stderr,"adaptive triangle delay to %d times the width\n",width_delay);
  }
  else
  {
    max_delay = delay;
    width_delay = 0;
    fprintf(stderr,"maximal triangle delay is %d\n",max_delay);
  }
  
  v_count = 0;
  f_count = 0;

  return true;
}

void SMwriter_smd::close(bool close_file, bool update_header)
{
  // close of SMwriter_smd
  while (triangle_buffer_size)
  {
    compress_triangle();
  }

  re_conn_op->encode(rmDone, 1); // done

  finishEncoder(v_count);
  finishModels();

  if (pq)
  {
    ic[0]->FinishCompressor();
    ic[1]->FinishCompressor();
    ic[2]->FinishCompressor();
  }
  else
  {
    fc[0]->FinishCompressor(0);
    fc[1]->FinishCompressor(0);
    fc[2]->FinishCompressor(0);
  }
  
  delete waiting_queue;
  delete traversal_queue;
  delete little_cache;

  if (dv->size()) fprintf(stderr,"WARNING: there are %d unfinalized vertices\n         this mesh will currently not decompress correctly\n",dv->size());
  if (vertex_hash->size()-dv->size()) fprintf(stderr,"WARNING: there are %d unused vertices\n         these vertices have not been compressed\n",vertex_hash->size()-dv->size());

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"none: %d last %d across %d\n", prediction_none, prediction_last, prediction_across);
  fprintf(stderr,"edge_buffer_size %d edge_buffer_maxsize %d\n",edge_buffer_size,edge_buffer_maxsize);
  fprintf(stderr,"vertex_buffer_size %d vertex_buffer_maxsize %d\n",vertex_buffer_size,vertex_buffer_maxsize);
  fprintf(stderr,"triangle_buffer_size %d triangle_buffer_maxsize %d\n",triangle_buffer_size,triangle_buffer_maxsize);
  fprintf(stderr,"right_confirm %d right_correct %d left_confirm %d left_correct %d\n",right_confirm,right_correct,left_confirm, left_correct);
  fprintf(stderr,"op_start %d op_add %d op_join %d op_fill %d op_end %d op_skip %d op_border %d\n",op_start,op_add,op_join,op_fill,op_end,op_skip,op_border);
  fprintf(stderr,"op_skip f_count *100 = %6.4f \n",100.0f*(float)op_skip/(float)f_count);

  fprintf(stderr,"%d %d max_in_width %d max_out_width %d diff %6.3f\n",vertex_hash->size(),dv->size(),max_in_width,max_out_width,100.0f*(max_out_width-max_in_width)/max_in_width);
  fprintf(stderr,"max_in_span %d max_out_span %d diff  %6.3f \n",max_in_span,max_out_span,100.0f*(max_out_span-max_in_span)/max_in_span);
#endif

  delete dv;
  delete vertex_hash;

  // close of SMwriter interface
  if (nverts != -1 && nverts != v_count)  fprintf(stderr,"WARNING: set nverts %d but v_count %d\n",nverts,v_count);
  if (nfaces != -1 && nfaces != f_count)  fprintf(stderr,"WARNING: set nfaces %d but f_count %d\n",nfaces,f_count);
  nverts = v_count;
  nfaces = f_count;
  v_count = -1;
  f_count = -1;
}

void SMwriter_smd::write_header()
{
  // write nverts
  if (nverts == -1)
  {
    re_conn->encode(2,0);
  }
  else
  {
    re_conn->encode(2,1);
    re_conn->encodeInt(nverts);
  }
  // write nfaces
  if (nfaces == -1)
  {
    re_conn->encode(2,0);
  }
  else
  {
    re_conn->encode(2,1);
    re_conn->encodeInt(nfaces);
  }
  // write bounding box
  if (bb_min_f == 0)
  {
    re_geom->encode(2,0);
  }
  else
  {
    re_geom->encode(2,1);
    re_geom->encodeFloat(bb_min_f[0]);
    re_geom->encodeFloat(bb_min_f[1]);
    re_geom->encodeFloat(bb_min_f[2]);
  }
  if (bb_max_f == 0)
  {
    re_geom->encode(2,0);
  }
  else
  {
    re_geom->encode(2,1);
    re_geom->encodeFloat(bb_max_f[0]);
    re_geom->encodeFloat(bb_max_f[1]);
    re_geom->encodeFloat(bb_max_f[2]);
  }
  // write comments
  if (true) // if have no comments
  {
    re_conn->encode(2,0);
  }
  else
  {
    /* yet to be implemented */
  }

  if (bb_min_f && bb_max_f && nbits)
  {
    if (true) // if want to use integer quantization
    {
      re_conn->encode(2,0);
    }
    else
    {
      /* yet to be implemented */
    }

    fprintf(stderr,"the bounding box is quantized to %d bits\n",nbits);

    pq = new PositionQuantizerNew();
    pq->SetMinMax(bb_min_f, bb_max_f);
    pq->SetPrecision(nbits);
    pq->SetupQuantizer();

    ic[0] = new IntegerCompressorNew();
    ic[0]->SetRange(pq->m_aiQuantRange[0]);
    ic[0]->SetPrecision(nbits);
    ic[0]->SetupCompressor(re_geom);

    ic[1] = new IntegerCompressorNew();
    ic[1]->SetRange(pq->m_aiQuantRange[1]);
    ic[1]->SetPrecision(nbits);
    ic[1]->SetupCompressor(re_geom);

    ic[2] = new IntegerCompressorNew();
    ic[2]->SetRange(pq->m_aiQuantRange[2]);
    ic[2]->SetPrecision(nbits);
    ic[2]->SetupCompressor(re_geom);
  }
  else
  {

    fprintf(stderr,"no bounding box! adaptive to %d bits\n",nbits);

    pq = 0;
    fc[0] = new FloatCompressor();
    fc[0]->SetPrecision(nbits);
    fc[0]->SetupCompressor(re_geom,0);

    fc[1] = new FloatCompressor();
    fc[1]->SetPrecision(nbits);
    fc[1]->SetupCompressor(re_geom,0);

    fc[2] = new FloatCompressor();
    fc[2]->SetPrecision(nbits);
    fc[2]->SetupCompressor(re_geom,0);
  }
}

SMwriter_smd::SMwriter_smd()
{
  // init of SMwriter_smd
}

SMwriter_smd::~SMwriter_smd()
{
  // clean-up for SMwriter_smd
}
