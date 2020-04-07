/*
===============================================================================

  FILE:  smwriter_smc.cpp
  
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
#include "smwriter_smc.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "rangemodel.h"
#include "rangeencoder.h"

#include "dynamicvector.h"
#include "littlecache.h"

#include "floatcompressor.h"

#include "positionquantizer_new.h"
#include "integercompressor_new.h"

#include "vec3fv.h"
#include "vec3iv.h"

#include "hash_map.h"

#define PRINT_CONTROL_OUTPUT
#undef PRINT_CONTROL_OUTPUT

#define SM_VERSION_SME 1
#define SM_VERSION_SME_NON_FINALIZED_EOF 3

#define SMC_START 0
#define SMC_ADD 1
#define SMC_JOIN 2
#define SMC_FILL 3
#define SMC_END 4
#define SMC_FILL_END 3

// data structures used for streaming compression

struct SMedge;

typedef struct SMvertex
{
  SMvertex* buffer_next;    // used for efficient memory management and by dynamicvector data structure
  float v[3];
  int index;
  int use_count;
  int list_size;
  int list_alloc;
  SMedge** list;
} SMvertex;

typedef struct SMedge
{
  SMedge* buffer_next;    // used for efficient memory management
  SMvertex* origin;
  SMvertex* target;
  float across[3];
} SMedge;

typedef hash_map<int, SMvertex*> my_vertex_hash;

static my_vertex_hash* vertex_hash;

static DynamicVector* dv;

static LittleCache* lc;

static FloatCompressor* fc[3];

static PositionQuantizerNew* pq;
static IntegerCompressorNew* ic[3];

// rangecoder and probability tables

#define MAX_DEGREE_ONE 4
#define MAX_USE_COUNT 15

static RangeEncoder* re_conn;
static RangeEncoder* re_conn_op;
static RangeEncoder* re_conn_cache;
static RangeEncoder* re_conn_index;
static RangeEncoder* re_conn_final;

static RangeEncoder* re_geom;

// is there more to encode
static RangeModel* rmDone;

// what was the last operation
static int last_op;

// codes next operation
static RangeModel** rmOp;

// codes non-manifoldness of start operations
static RangeModel* rmS_Old;

// codes cache hits for start operations
static RangeModel** rmS_Cache;

// codes cache hits for add/join operations
static RangeModel** rmAJ_Cache;

// codes cache hits fill/end operations
static RangeModel** rmFE_Cache;

// codes vertex finalization
static RangeModel*** rmFinalized;

// statistics

#ifdef PRINT_CONTROL_OUTPUT
static int op_start = 0;
static int op_add = 0;
static int op_join = 0;
static int op_fill_end = 0;
static int op_fill = 0;
static int op_end = 0;

static int used_index = 0;
static int used_cache = 0;

static int add_miss = 0;
static int add_hit[] = {0,0,0,0,0,0};

static int fill_miss = 0;
static int fill_hit[] = {0,0,0,0,0,0,0,0,0};

static int prediction_none = 0;
static int prediction_last = 0;
static int prediction_across = 0;

static int vertex_buffer_maxsize = 0;
static int edge_buffer_maxsize = 0;
#endif

static void initEncoder(FILE* file)
{
  if (file)
  {
    re_conn = new RangeEncoder(file);
    re_conn_op = re_conn;
    re_conn_cache = re_conn;
    re_conn_index = re_conn;
    re_conn_final = re_conn;
    re_geom = re_conn;
  }
  else
  {
    re_conn = new RangeEncoder(0,false);
    re_conn_op = new RangeEncoder(0,false);
    re_conn_cache = new RangeEncoder(0,false);
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
    re_conn_cache->done();
    re_conn_index->done();
    re_conn_final->done();
    re_conn->done();
    re_geom->done();

    fprintf(stderr,"\n");
    fprintf(stderr,"*** compression rate details ***\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"op:   \t%6.3f bpv\n", 8.0f/nverts*re_conn_op->getNumberBytes());
    fprintf(stderr,"cache:\t%6.3f bpv\n", 8.0f/nverts*re_conn_cache->getNumberBytes());
    fprintf(stderr,"index:\t%6.3f bpv\n", 8.0f/nverts*re_conn_index->getNumberBytes());
    fprintf(stderr,"final:\t%6.3f bpv\n", 8.0f/nverts*re_conn_final->getNumberBytes());
    fprintf(stderr,"rest: \t%6.3f bpv\n", 8.0f/nverts*re_conn->getNumberBytes());
    fprintf(stderr,"\n");
    fprintf(stderr,"conn: \t%6.3f bpv\n", 8.0f/nverts*(re_conn_op->getNumberBytes()+re_conn_cache->getNumberBytes()+re_conn_index->getNumberBytes()+re_conn_final->getNumberBytes()+re_conn->getNumberBytes()));
    fprintf(stderr,"geom: \t%6.3f bpv\n", 8.0f/nverts*re_geom->getNumberBytes());
    fprintf(stderr,"\n");
    fprintf(stderr,"total:\t%6.3f bpv\n", 8.0f/nverts*(re_geom->getNumberBytes()+re_conn_op->getNumberBytes()+re_conn_cache->getNumberBytes()+re_conn_index->getNumberBytes()+re_conn_final->getNumberBytes()+re_conn->getNumberBytes()));

#ifdef PRINT_CONTROL_OUTPUT
    fprintf(stderr,"none: %d last %d across %d\n", prediction_none, prediction_last, prediction_across);
    if (pq)
    {
      fprintf(stderr,"small %d %d %d %f\n", ic[0]->num_predictions_small, ic[1]->num_predictions_small, ic[2]->num_predictions_small, 100.0f*(ic[0]->num_predictions_small+ic[1]->num_predictions_small+ic[2]->num_predictions_small)/3/nverts);
    }
#endif
    
    delete re_conn;
    delete re_conn_op;
    delete re_conn_cache;
    delete re_conn_index;
    delete re_conn_final;
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
  rmOp = (RangeModel**)malloc(sizeof(RangeModel*)*5);
  rmOp[0] = new RangeModel(4,0,compress); // after SMC_START
  rmOp[1] = new RangeModel(4,0,compress); // after SMC_ADD
  rmOp[2] = new RangeModel(4,0,compress); // after SMC_JOIN
  rmOp[3] = new RangeModel(4,0,compress); // after SMC_FILL
  rmOp[4] = new RangeModel(4,0,compress); // after SMC_END

  // range tables for start (S) operation
  rmS_Old = new RangeModel(5,0,compress); // number of old (e.g. non-manifold) vertices in a start operation (4 = NON_FINALIZED_EOF)
  rmS_Cache = (RangeModel**)malloc(sizeof(RangeModel*)*3);
  rmS_Cache[0] = new RangeModel(4,0,compress); // for 0
  rmS_Cache[1] = new RangeModel(4,0,compress); // for 1
  rmS_Cache[2] = new RangeModel(4,0,compress); // for 2

  // range tables for add (A) / join (J) operations

  rmAJ_Cache = (RangeModel**)malloc(sizeof(RangeModel*)*5);
  rmAJ_Cache[0] = new RangeModel(7,0,compress); // after SMC_START
  rmAJ_Cache[1] = new RangeModel(7,0,compress); // after SMC_ADD
  rmAJ_Cache[2] = new RangeModel(7,0,compress); // after SMC_JOIN
  rmAJ_Cache[3] = new RangeModel(7,0,compress); // after SMC_FILL
  rmAJ_Cache[4] = new RangeModel(7,0,compress); // after SMC_END

  // range tables for fill (F) / end (E) operations

  rmFE_Cache = (RangeModel**)malloc(sizeof(RangeModel*)*5);
  rmFE_Cache[0] = new RangeModel(10,0,compress); // after SMC_START
  rmFE_Cache[1] = new RangeModel(10,0,compress); // after SMC_ADD
  rmFE_Cache[2] = new RangeModel(10,0,compress); // after SMC_JOIN
  rmFE_Cache[3] = new RangeModel(10,0,compress); // after SMC_FILL
  rmFE_Cache[4] = new RangeModel(10,0,compress); // after SMC_END

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

  // range table for whether there is more
  delete rmDone;

  // range table for which operation it is
  delete rmOp[0];
  delete rmOp[1];
  delete rmOp[2];
  delete rmOp[3];
  delete rmOp[4];
  free(rmOp);

  // range tables for start (S) operation
  delete rmS_Old;

  // range tables for start (S) operations
  delete rmS_Cache[0];
  delete rmS_Cache[1];
  delete rmS_Cache[2];
  free(rmS_Cache);

  // range tables for add (A) / join (J) operations
  delete rmAJ_Cache[0];
  delete rmAJ_Cache[1];
  delete rmAJ_Cache[2];
  delete rmAJ_Cache[3];
  delete rmAJ_Cache[4];
  free(rmAJ_Cache);

  // range tables for fill (F) / end (E) operations
  delete rmFE_Cache[0];
  delete rmFE_Cache[1];
  delete rmFE_Cache[2];
  delete rmFE_Cache[3];
  delete rmFE_Cache[4];
  free(rmFE_Cache);

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

// efficient memory allocation

static int vertex_buffer_size = 0;
static int vertex_buffer_alloc = 16;
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
      vertex_buffer_next[i].list = 0;
    }
    vertex_buffer_next[vertex_buffer_alloc-1].buffer_next = 0;
    vertex_buffer_alloc = 2*vertex_buffer_alloc;
  }
  // get pointer to next available vertex
  SMvertex* vertex = vertex_buffer_next;
  vertex_buffer_next = vertex->buffer_next;

  if (vertex->list == 0)
  {
    vertex->list = (SMedge**)malloc(sizeof(SMedge*)*10);
    vertex->list_alloc = 10;
  }
  vertex->list_size = 0;
  vertex->use_count = 0;

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
static int edge_buffer_alloc = 16;
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

static void compressVertexPosition(const float* l, float* n)
{
#ifdef PRINT_CONTROL_OUTPUT
  prediction_last++;
#endif

  if (pq)
  {
    pq->EnQuantize(n, (int*)n);
    for (int i = 0; i < 3; i++)
    {
      ic[i]->CompressLast(((const int*)l)[i],((int*)n)[i]);
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

static void addEdgeToVertices(SMedge* edge)
{
  addEdgeToVertex(edge, edge->origin);
  addEdgeToVertex(edge, edge->target);
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

void SMwriter_smc::write_vertex(const float* v_pos_f)
{
  if (v_count == 0) write_header();

  SMvertex* vertex = allocVertex();
  vertex->index = v_count;
  VecCopy3fv(vertex->v, v_pos_f);
  vertex_hash->insert(my_vertex_hash::value_type(v_count, vertex));
  v_count++;
}

void SMwriter_smc::write_triangle(const int* t_idx)
{
  fprintf(stderr, "ERROR: write_triangle(const int* t_idx) not supported by SMwriter_smc\n");
  exit(0);
}

void SMwriter_smc::write_triangle(const int* t_idx, const bool* t_final)
{
  int i,j;
  int dv_index;
  int lc_pos;
  int candidate;
  int candidate_count;
  int degree_one;
  int use_count;
  int rot = 0;

  my_vertex_hash::iterator hash_elements[3];
  SMvertex* vertices[3];
  SMedge* edges[3];

  // get vertices from hash and check if vertices already visited

  SMvertex* v_visited[3];
  int num_v_visited = 0;
  int last_v_visited;
  int last_v_unvisited;
  for (i = 0; i < 3; i++)
  {
    hash_elements[i] = vertex_hash->find(t_idx[i]);
    if (hash_elements[i] == vertex_hash->end())
    {
      fprintf(stderr,"ERROR: vertex %d not in hash\n",t_idx[i]);
      exit(0);
    }

    v_visited[i] = (*hash_elements[i]).second;

    if (v_visited[i]->use_count)
    {
      last_v_visited = i;
      num_v_visited++;
    }
    else
    {
      v_visited[i] = 0;
      last_v_unvisited = i;
    }
  }

  // are edges already visited once before with opposite orientation

  SMedge* e_visited[3];
  int num_e_visited = 0;
  int last_e_visited = -1;
  int last_e_unvisited;

  for (i = 0; i < 3; i++)
  {
    e_visited[i] = 0;
    if (v_visited[i] && v_visited[(i+1)%3])
    {
      for (j = 0; j < v_visited[i]->list_size; j++)
      {
        if (v_visited[i]->list[j]->origin == v_visited[(i+1)%3])
        {
          e_visited[i] = v_visited[i]->list[j];
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

  // compress the triangle based in its operation: start, add, join, or fill/end

  if (num_e_visited == 0) // start
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_start++;
#endif

    if (edge_buffer_size)
    {
      re_conn_op->encode(rmOp[last_op], SMC_START);
    }
    else
    {
      re_conn_op->encode(rmDone, 0); // more to come
    }
    last_op = SMC_START;

    // specify how many old (non-manifold) vertices there are

    re_conn_op->encode(rmS_Old, num_v_visited);

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

    vertices[0] = (*hash_elements[rot]).second;
    vertices[1] = (*hash_elements[(rot+1)%3]).second;
    vertices[2] = (*hash_elements[(rot+2)%3]).second;

    edges[0] = 0;
    edges[1] = 0;
    edges[2] = 0;

    // encode vertex v0

    if (vertices[0]->use_count)
    {
      // encode the index of v0 (-> improve with cache of with prediction)
      lc_pos = lc->pos(vertices[0]->index);
      if (lc_pos == -1)
      {
        re_conn_cache->encode(rmS_Cache[0],3);        
        // encode its index among all active vertices
        dv_index = dv->getRelativeIndex(vertices[0]);
        re_conn_index->encode(dv->size(),dv_index);        
#ifdef PRINT_CONTROL_OUTPUT
        used_index++;
#endif
      }
      else
      {
        re_conn_cache->encode(rmS_Cache[0],lc_pos);        
#ifdef PRINT_CONTROL_OUTPUT
        used_cache++;
#endif
      }
    }
    else
    {
      // encode its position
      compressVertexPosition(vertices[0]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[0]);
    }

    // encode vertex v1

    if (vertices[1]->use_count)
    {
      // encode the index of v1 (-> improve with cache of with prediction)
      lc_pos = lc->pos(vertices[1]->index);
      if (lc_pos == -1)
      {
        re_conn_cache->encode(rmS_Cache[1],3);        
        // encode its index among all active vertices
        dv_index = dv->getRelativeIndex(vertices[1]);
        re_conn_index->encode(dv->size(),dv_index);        
#ifdef PRINT_CONTROL_OUTPUT
        used_index++;
#endif
      }
      else
      {
        re_conn_cache->encode(rmS_Cache[1],lc_pos);        
#ifdef PRINT_CONTROL_OUTPUT
        used_cache++;
#endif
      }
    }
    else
    {
      // encode its position
      compressVertexPosition(vertices[0]->v, vertices[1]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[1]);
    }

    // encode vertex v2

    if (vertices[2]->use_count)
    {
      // encode the index of v2 (-> improve with cache of with prediction)
      lc_pos = lc->pos(vertices[2]->index);
      if (lc_pos == -1)
      {
        re_conn_cache->encode(rmS_Cache[2],3);        
        // encode its index among all active vertices
        dv_index = dv->getRelativeIndex(vertices[2]);
        re_conn_index->encode(dv->size(),dv_index);        
#ifdef PRINT_CONTROL_OUTPUT
        used_index++;
#endif
      }
      else
      {
        re_conn_cache->encode(rmS_Cache[2],lc_pos);        
#ifdef PRINT_CONTROL_OUTPUT
        used_cache++;
#endif
      }
    }
    else
    {
      // encode its position
      compressVertexPosition(vertices[0]->v, vertices[2]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[2]);
    }
  }
  else if (num_e_visited == 1) // add or join)
  {
    if (num_v_visited == 2) // add
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_add++;
#endif
      re_conn_op->encode(rmOp[last_op], SMC_ADD);
    }
    else
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_join++;
#endif
      re_conn_op->encode(rmOp[last_op], SMC_JOIN);
    }

    // rotate triangle if necessary

    rot = last_e_visited;

    // encode the (rotated) add (or join) configuration

    vertices[0] = v_visited[rot];
    vertices[1] = v_visited[(rot+1)%3];
    vertices[2] = (*hash_elements[(rot+2)%3]).second;

    edges[0] = e_visited[rot];
    edges[1] = 0;
    edges[2] = 0;

    // try to encode the index of v0 with cache
    lc_pos = lc->pos(vertices[0]->index);
    if (lc_pos == -1)
    {
      // try to encode the index of v1 with cache
      lc_pos = lc->pos(vertices[1]->index);
      if (lc_pos == -1)
      {
        re_conn_cache->encode(rmAJ_Cache[last_op],6);        
        dv_index = dv->getRelativeIndex(vertices[0]);
        re_conn_index->encode(dv->size(),dv_index);        
#ifdef PRINT_CONTROL_OUTPUT
        add_miss++;
        used_index++;
#endif
      }
      else
      {
        lc_pos += 3;
        re_conn_cache->encode(rmAJ_Cache[last_op],lc_pos);
#ifdef PRINT_CONTROL_OUTPUT
        add_hit[lc_pos]++;
        used_cache++;
#endif
      }
    }
    else
    {
      re_conn_cache->encode(rmAJ_Cache[last_op],lc_pos);        
#ifdef PRINT_CONTROL_OUTPUT
      add_hit[lc_pos]++;
      used_cache++;
#endif
    }

    if (lc_pos < 3)
    {
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
      }
    }
    else
    {
      // encode which of the edges incident to v1 is e0
      candidate_count = 0;

      // how many potential candidates are around that vertex?
      for (i = 0; i < vertices[1]->list_size; i++)
      {
        if (vertices[1]->list[i]->origin == vertices[1])
        {
          if (vertices[1]->list[i] == edges[0]) candidate = candidate_count;
          candidate_count++;
        }
      }
      // only if there is more than one we need to encode it
      if (candidate_count > 1)
      {
        re_conn->encode(candidate_count,candidate);
      }
    }

    // encode v2
    if (vertices[2]->use_count) // join
    {
      // encode its index among all active vertices
      dv_index = dv->getRelativeIndex(vertices[2]);
      re_conn_index->encode(dv->size(),dv_index);        
#ifdef PRINT_CONTROL_OUTPUT
      used_index++;
#endif
      last_op = SMC_JOIN;
    }
    else // add
    {
      // encode its position
      compressVertexPosition(vertices[0]->v, edges[0]->across, vertices[1]->v, vertices[2]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[2]);
      last_op = SMC_ADD;
    }
  }
  else // fill or end
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_fill_end++;
#endif
    re_conn_op->encode(rmOp[last_op], SMC_FILL_END);

    // rotate triangle if necessary

    if (num_e_visited == 2)
    {
      rot = (last_e_unvisited + 2) % 3;
    }
    else
    {
      for (i = 0; i < 3; i++) // try to find a vertex that is in the cache
      {
        if (lc->pos(v_visited[i]->index) != -1)
        {
          rot = i;
          break;
        }
      }
    }

    // encode the (rotated) fill (or end) configuration

    vertices[0] = v_visited[rot];
    vertices[1] = v_visited[(rot+1)%3];
    vertices[2] = v_visited[(rot+2)%3];

    edges[0] = e_visited[rot];
    edges[1] = e_visited[(rot+1)%3];
    edges[2] = e_visited[(rot+2)%3];

    // see whether either v0, v1, or v2 is a cached vertex

    lc_pos = lc->pos(vertices[0]->index);
    if (lc_pos == -1)
    {
      lc_pos = lc->pos(vertices[1]->index);
      if (lc_pos == -1)
      {
        lc_pos = lc->pos(vertices[2]->index);
        if (lc_pos == -1)
        {
          // neither vertex is in the cache ... using vertex v0
          re_conn_cache->encode(rmFE_Cache[last_op],9);
          dv_index = dv->getRelativeIndex(vertices[0]);
          re_conn_index->encode(dv->size(),dv_index);
#ifdef PRINT_CONTROL_OUTPUT
          fill_miss++;
          used_index++;
#endif
        }
        else
        {
          lc_pos += 6;
          re_conn_cache->encode(rmFE_Cache[last_op],lc_pos);
#ifdef PRINT_CONTROL_OUTPUT
          fill_hit[lc_pos]++;
          used_cache++;
#endif
        }
      }
      else
      {
        lc_pos += 3;
        re_conn_cache->encode(rmFE_Cache[last_op],lc_pos);        
#ifdef PRINT_CONTROL_OUTPUT
        fill_hit[lc_pos]++;
        used_cache++;
#endif
      }
    }
    else
    {
      re_conn_cache->encode(rmFE_Cache[last_op],lc_pos);        
#ifdef PRINT_CONTROL_OUTPUT
      fill_hit[lc_pos]++;
      used_cache++;
#endif
    }

    if (lc_pos < 3) // we have v0
    {
      // encode edge e0 (e.g. the edge from v0 to v1)
      candidate_count = 0;

      // how many edges are around v0 that could potentially be e0?
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
      }

      // encode edge e2 (e.g. the edge from v2 to v0)
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
      }
    }
    else if (lc_pos < 6) // we have v1
    {
      // encode edge e0 (e.g. the edge from v0 to v1)
      candidate_count = 0;

      // how many edges are around v1 that could potentially be e0?
      for (i = 0; i < vertices[1]->list_size; i++)
      {
        if (vertices[1]->list[i]->origin == vertices[1])
        {
          if (vertices[1]->list[i] == edges[0]) candidate = candidate_count;
          candidate_count++;
        }
      }

      // only if there is more than one we need to encode it
      if (candidate_count > 1)
      {
        re_conn->encode(candidate_count,candidate);
      }

      // encode edge e2 (e.g. the edge from v2 to v0)
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
      }
    }
    else // we have v2
    {
      // encode edge e2 (e.g. the edge from v2 to v0)
      candidate_count = 0;

      // how many edges are around v2 that could potentially be e2?
      for (i = 0; i < vertices[2]->list_size; i++)
      {
        if (vertices[2]->list[i]->target == vertices[2])
        {
          if (vertices[2]->list[i] == edges[2]) candidate = candidate_count;
          candidate_count++;
        }
      }

      // only if there is more than one we need to encode it
      if (candidate_count > 1)
      {
        re_conn->encode(candidate_count,candidate);
      }

      // encode edge e0 (e.g. the edge from v0 to v1)
      candidate_count = 0;

      // how many edges are around v0 that could potentially be e0?
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
      }
    }

    if (num_e_visited == 2)
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_fill++;
#endif
      last_op = SMC_FILL;
    }
    else
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_end++;
#endif
      last_op = SMC_END;

      // make sure that the decoder uses the same edges[1]
      for (i = 0; i < vertices[2]->list_size; i++)
      {
        if (vertices[2]->list[i]->target == vertices[1])
        {
          edges[1] = vertices[2]->list[i];
          break;
        }
      }
    }
  }

  // update cache

  lc->put(vertices[0]->index, vertices[1]->index, vertices[2]->index, vertices[0], vertices[1], vertices[2]);

  // increment vertex use_counts, create edges, and update edge degrees

  for (i = 0; i < 3; i++)
  {
    vertices[i]->use_count++;

    if (edges[i] == 0)
    {
      edges[i] = allocEdge(vertices[(i+2)%3]->v);
      edges[i]->origin = vertices[i];
      edges[i]->target = vertices[(i+1)%3];
      addEdgeToVertices(edges[i]);
    }
    else
    {
      removeEdgeFromVertices(edges[i]);
      deallocEdge(edges[i]);
    }
  }

  // write finalization info and finalize vertices

  for (i = 0; i < 3; i++)
  {
    j = (i+rot)%3;

    degree_one = vertices[i]->list_size;
    if (degree_one >= MAX_DEGREE_ONE)
    {
      degree_one = MAX_DEGREE_ONE-1;
    }
    use_count = vertices[i]->use_count;
    if (use_count >= MAX_USE_COUNT)
    {
      use_count = MAX_USE_COUNT-1;
    }
    re_conn_final->encode(rmFinalized[degree_one][use_count],t_final[j]);

    if (t_final[j])
    {
      SMvertex* vertex = vertices[i];
      vertex_hash->erase(hash_elements[j]);
      dv->removeElement(vertex);
      for (j = 0; j < vertex->list_size; j++)
      {
        SMedge* edge = vertex->list[j];
        if (vertex == edge->origin)
        {
          removeEdgeFromVertex(edge, edge->target);
        }
        else
        {
          removeEdgeFromVertex(edge, edge->origin);
        }
        deallocEdge(edge);
      }
      deallocVertex(vertex);
    }
  }
  f_count++;
}

void SMwriter_smc::write_finalized(int final_idx)
{
  fprintf(stderr, "ERROR: write_finalized(int final_idx) not supported by SMwriter_smc\n");
  exit(0);
}

bool SMwriter_smc::open(FILE* file, int nbits)
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
    fputc(SM_VERSION_SME_NON_FINALIZED_EOF, file);
  }

  initVertexBuffer(1024);
  initEdgeBuffer(1024);

  vertex_hash = new my_vertex_hash;
  dv = new DynamicVector();
  lc = new LittleCache();

  initEncoder(file);
  initModels(1);

  last_op = 0;

  // write precision
  re_conn->encode(25,nbits);

  // store precision 
  this->nbits = nbits;

  v_count = 0;
  f_count = 0;

  return true;
}

void SMwriter_smc::close(bool close_file, bool update_header)
{
  // close of SMwriter_smc
  if (edge_buffer_size == 0)
  {
    re_conn_op->encode(rmDone, 1); // done
  }
  else
  {
    re_conn_op->encode(rmOp[last_op], SMC_START);
    re_conn_op->encode(rmS_Old, 4);
  }
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

  if (dv->size()) fprintf(stderr,"WARNING: there are %d unfinalized vertices\n",dv->size());
  if (vertex_hash->size()-dv->size()) fprintf(stderr,"WARNING: %d unused vertices have not been compressed\n",vertex_hash->size()-dv->size());

  delete dv;
  delete vertex_hash;

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"nfaces %d f_count %d ops %d\n",nfaces,f_count,op_start+op_add+op_join+op_fill_end);
  fprintf(stderr,"edge_buffer_size %d edge_buffer_maxsize %d\n",edge_buffer_size,edge_buffer_maxsize);
  fprintf(stderr,"vertex_buffer_size %d vertex_buffer_maxsize %d\n",vertex_buffer_size,vertex_buffer_maxsize);
  fprintf(stderr,"op_start %d (%4.2f) op_add %d (%4.2f) od_join %d (%4.2f) op_fill %d (%4.2f) op_end %d (%4.2f)\n",op_start,100.0f*op_start/(op_start+op_add+op_join+op_fill_end),op_add,100.0f*op_add/(op_start+op_add+op_join+op_fill_end),op_join,100.0f*op_join/(op_start+op_add+op_join+op_fill_end),op_fill,100.0f*op_fill/(op_start+op_add+op_join+op_fill_end),op_end,100.0f*op_end/(op_start+op_add+op_join+op_fill_end));
  fprintf(stderr,"add_miss %d add_hit %d (%d %d %d %d %d %d)\n",add_miss,add_hit[0]+add_hit[1]+add_hit[2]+add_hit[3]+add_hit[4]+add_hit[5],add_hit[0],add_hit[1],add_hit[2],add_hit[3],add_hit[4],add_hit[5]);
  fprintf(stderr,"fill_miss %d fill_hit %d (%d %d %d %d %d %d %d %d %d)\n",fill_miss,fill_hit[0]+fill_hit[1]+fill_hit[2]+fill_hit[3]+fill_hit[4]+fill_hit[5]+fill_hit[6]+fill_hit[7]+fill_hit[8],fill_hit[0],fill_hit[1],fill_hit[2],fill_hit[3],fill_hit[4],fill_hit[5],fill_hit[6],fill_hit[7],fill_hit[8]);
  fprintf(stderr,"used_index %d (%4.2f) used_cache %d  (%4.2f)\n",used_index,100.0f*used_index/(used_index+used_cache),used_cache,100.0f*used_cache/(used_index+used_cache));
#endif

  // close of SMwriter interface
  if (nverts != -1 && nverts != v_count)  fprintf(stderr,"WARNING: set nverts %d but v_count %d\n",nverts,v_count);
  if (nfaces != -1 && nfaces != f_count)  fprintf(stderr,"WARNING: set nfaces %d but f_count %d\n",nfaces,f_count);
  nverts = v_count;
  nfaces = f_count;
  v_count = -1;
  f_count = -1;
}

void SMwriter_smc::write_header()
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

  // setup geometry compressor

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

#ifdef PRINT_CONTROL_OUTPUT
    fprintf(stderr,"the bounding box is quantized to %d bits\n",nbits);
#endif

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

#ifdef PRINT_CONTROL_OUTPUT
    fprintf(stderr,"no bounding box! adaptive to %d bits\n",nbits);
#endif

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

SMwriter_smc::SMwriter_smc()
{
  // re-init of SMwriter interface
  need_pre_order = true;
  // init of SMwriter_smc
}

SMwriter_smc::~SMwriter_smc()
{
  // clean-up for SMwriter_smc
}
