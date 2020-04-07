/*
===============================================================================

  FILE:  smwriter_smc_old.cpp
  
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
#include "smwriter_smc_old.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "littlecache.h"
#include "dynamicvector.h"

#include "rangemodel.h"
#include "rangeencoder.h"

#include "floatcompressor.h"

#include "positionquantizer.h"
#include "integercompressor.h"

#include "vec3fv.h"
#include "vec3iv.h"

#include "hash_map.h"

//#define PRINT_CONTROL_OUTPUT fprintf
#define PRINT_DEBUG_OUTPUT if (0) fprintf
#define PRINT_CONTROL_OUTPUT if (0) fprintf

#define USE_CACHE 1

#define SMC_ADD_JOIN 0
#define SMC_FILL_END 1
#define SMC_START 2

// data structures used for streaming compression

struct SMedge;

typedef struct SMvertex
{
  int dynamicvector; // used by dynamicvector data structure
  SMvertex* buffer_next;    // used for efficient memory management
  float v[3];
  int index;
  int use_count;
  int degree_one;
  int list_size;
  int list_alloc;
  SMedge** list;
} SMvertex;

typedef struct SMedge
{
  SMedge* buffer_next;    // used for efficient memory management
  SMvertex* origin;
  SMvertex* target;
  int degree;
  float across[3];
} SMedge;

typedef hash_map<int, SMvertex*> my_vertex_hash;

static my_vertex_hash* vertex_hash;

static DynamicVector* dv;

static LittleCache* lc;
 
static int precision_bits;

static FloatCompressor* fc[3];

static PositionQuantizer* pq;
static IntegerCompressor* ic[3];

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

// which operation
static int last_SAJFE = 0;
static RangeModel** rmSAJFE;

// used for start operations
static RangeModel* rmS_0Old;
static RangeModel* rmS_1Old;
static RangeModel* rmS_2Old;

//RangeModel* rmS_CacheHit;
//RangeModel* rmS_CachePos;

// used for add/join operations
static RangeModel* rmAJ_Manifold;
static RangeModel* rmAJ_ManifoldOriented;
static RangeModel* rmAJ_NonManifoldOriented;
static RangeModel* rmAJ_Old;

static RangeModel* rmAJ_CacheHit;
static RangeModel* rmAJ_CachePos;

// used for fill/end operations
static RangeModel* rmFE_0Manifold;
static RangeModel* rmFE_0ManifoldOriented;
static RangeModel* rmFE_0NonManifoldOriented;
static RangeModel* rmFE_2Manifold;
static RangeModel* rmFE_2ManifoldOriented;
static RangeModel* rmFE_2NonManifoldOriented;

static RangeModel* rmFE_CacheHit;
static RangeModel* rmFE_CachePos;

// used for vertex finalization
static RangeModel*** rmFinalized;

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
    PRINT_CONTROL_OUTPUT(stderr,"total: bytes %d bpv %f\n", re_conn->getNumberBytes(),(float)re_conn->getNumberBits()/nverts);
    delete re_conn;
  }
}

static void initModels(int compress)
{
  int i,j;

  // range table for whether there is more
  rmDone = new RangeModel(2,0,compress);

  // range table for which operation it is
  rmSAJFE = (RangeModel**)malloc(sizeof(RangeModel*)*3);
  rmSAJFE[0] = new RangeModel(3,0,compress);
  rmSAJFE[1] = new RangeModel(3,0,compress);
  rmSAJFE[2] = new RangeModel(3,0,compress);

  // range tables for start (S) operation
  rmS_0Old = new RangeModel(2,0,compress);
  rmS_1Old = new RangeModel(2,0,compress);
  rmS_2Old = new RangeModel(2,0,compress);

//  rmS_CacheHit = new RangeModel(2,0,compress);
//  rmS_CachePos = new RangeModel(3,0,compress);

  // range tables for add (A) / join (J) operations
  rmAJ_Manifold = new RangeModel(2,0,compress);
  rmAJ_ManifoldOriented = new RangeModel(2,0,compress);
  rmAJ_NonManifoldOriented = new RangeModel(2,0,compress);
  rmAJ_Old = new RangeModel(2,0,compress);

  rmAJ_CacheHit = new RangeModel(2,0,compress);
  rmAJ_CachePos = new RangeModel(3,0,compress);

  // range tables for fill (F) / end (E) operations
  rmFE_0Manifold = new RangeModel(2,0,compress);
  rmFE_0ManifoldOriented = new RangeModel(2,0,compress);
  rmFE_0NonManifoldOriented = new RangeModel(2,0,compress);
  rmFE_2Manifold = new RangeModel(2,0,compress);
  rmFE_2ManifoldOriented = new RangeModel(2,0,compress);
  rmFE_2NonManifoldOriented = new RangeModel(2,0,compress);

  rmFE_CacheHit = new RangeModel(2,0,compress);
  rmFE_CachePos = new RangeModel(3,0,compress);

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
}

// statistics

static int op_start = 0;
static int op_add_join = 0;
static int op_fill_end = 0;

static int add_hit = 0;
static int add_miss = 0;
static int add_miss_hit_possible = 0;
static int add_pos[] = {0,0,0};

static int fill_hit = 0;
static int fill_miss = 0;
static int fill_miss_hit_possible_v1 = 0;
static int fill_miss_hit_possible_v2 = 0;
static int fill_pos[] = {0,0,0};

static int start_non_manifold = 0;
static int add_non_manifold = 0;

//static int output_type0;
//static int output_type1;
//static int output_type2;
//static int output_type3;

static int vertex_buffer_size;
static int edge_buffer_size;

static int vertex_buffer_maxsize;
static int edge_buffer_maxsize;

// efficient memory allocation

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
  vertex->degree_one = 0;

  vertex_buffer_size++; if (vertex_buffer_size > vertex_buffer_maxsize) vertex_buffer_maxsize = vertex_buffer_size;

  return vertex;
}

static void deallocVertex(SMvertex* vertex)
{
  vertex->buffer_next = vertex_buffer_next;
  vertex_buffer_next = vertex;
  vertex_buffer_size--;
}

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
  edge->degree = 1;
  VecCopy3fv(edge->across,v);

  edge_buffer_size++; if (edge_buffer_size > edge_buffer_maxsize) edge_buffer_maxsize = edge_buffer_size;

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
  if (pq)
  {
    int* qn = (int*)n;
    pq->EnQuantize(n, qn);
    for (int i = 0; i < 3; i++)
    {
      ic[i]->CompressNone(qn[i]);
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
  if (pq)
  {
    int* ql = (int*)l;
    int* qn = (int*)n;
    pq->EnQuantize(n, qn);
    int corr[3];
    pq->GetCorrector(ql, qn, corr); 
    for (int i = 0; i < 3; i++)
    {
      ic[i]->CompressLast(corr[i]);
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

static void compressVertexPosition(float* a, float* b, float* c, float* n)
{
  if (pq)
  {
    int* qa = (int*)a;
    int* qb = (int*)b;
    int* qc = (int*)c;
    int* qn = (int*)n;
    pq->EnQuantize(n, qn);
    int pred[3];
    VecAdd3iv(pred, qa, qc);
    VecSelfSubtract3iv(pred, qb);
    int corr[3];
    pq->GetCorrector(pred, qn, corr); 
    for (int i = 0; i < 3; i++)
    {
      ic[i]->CompressAcross(corr[i]);
    }
  }
  else
  {
    float pp[3];
    VecAdd3fv(pp, a, c);
    VecSelfSubtract3fv(pp, b);
    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->CompressAcross(pp[i],n[i]);
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
  vertex->degree_one++;
  vertex->list_size++;
}

static void removeEdgeFromVertex(SMedge* edge, SMvertex* vertex)
{
  int i,j;
  for (i = 0; i < vertex->list_size; i++)
  {
    if (vertex->list[i] == edge)
    {
      for (j = i+1; j < vertex->list_size; j++,i++)
      {
        vertex->list[i] = vertex->list[j];
      }
      vertex->list_size--;
      if (edge->degree == 1)
      {
        vertex->degree_one--;
      }
      return;
    }
  }
  fprintf(stderr,"ERROR: did not find edge in vertex list\n");
  exit(0);
}

void SMwriter_smc_old::write_vertex(const float* v_pos_f)
{
  if (v_count == 0) write_header();

  SMvertex* vertex = allocVertex();
  vertex->index = v_count;
  VecCopy3fv(vertex->v, v_pos_f);
  vertex_hash->insert(my_vertex_hash::value_type(v_count, vertex));
  v_count++;
}

void SMwriter_smc_old::write_triangle(const int* t_idx)
{
  fprintf(stderr, "ERROR: write_triangle(const int* t_idx) not supported by SMwriter_smc_old\n");
  exit(0);
}

void SMwriter_smc_old::write_triangle(const int* t_idx, const bool* t_final)
{
  int i,j;
  my_vertex_hash::iterator hash_elements[3];
  SMvertex* vertices[3];
  SMedge* edges[3];

  // get vertices from hash
  for (i = 0; i < 3; i++)
  {
    hash_elements[i] = vertex_hash->find(t_idx[i]);
    if (hash_elements[i] == vertex_hash->end())
    {
      fprintf(stderr,"ERROR: vertex %d not in hash\n",t_idx[i]);
      exit(0);
    }
    else
    {
      vertices[i] = (*hash_elements[i]).second;
    }
  }

  // are vertices already visited
  int v_visited[3];
  int num_v_visited = 0;
  int last_v_visited = -1;
  int last_v_unvisited = -1;

  for (i = 0; i < 3; i++)
  {
    if (vertices[i]->use_count)
    {
      v_visited[i] = 1;
      last_v_visited = i;
      num_v_visited++;
    }
    else
    {
      v_visited[i] = 0;
      last_v_unvisited = i;
    }
  }

  // are edges already visited
  int e_visited[3];
  int num_e_visited = 0;
  int last_e_visited = -1;
  int last_e_unvisited = -1;

  for (i = 0; i < 3; i++)
  {
    e_visited[i] = 0;
    if (v_visited[i] && v_visited[(i+1)%3])
    {
      for (j = 0; j < vertices[i]->list_size; j++)
      {
        if (vertices[i]->list[j]->origin == vertices[(i+1)%3])
        {
          e_visited[i] = 1;
          last_e_visited = i;
          edges[i] = vertices[i]->list[j];
          break;
        }
        else if (vertices[i]->list[j]->target == vertices[(i+1)%3]) // not-oriented case
        {
          e_visited[i] = 1;
          last_e_visited = i;
          edges[i] = vertices[i]->list[j];
          break;
        }
      }
      if (last_e_visited != i)
      {
        last_e_unvisited = i;
      }
      else
      {
        num_e_visited++;
      }
    }
  }

  // compress the triangle based in its operation: start, add/join, or fill/end

  int dv_index;
  int lc_pos;
  int candidate;
  int candidate_count;
  int degree_one;
  int use_count;
  int rot = 0;

  SMvertex* v0;
  SMvertex* v1;
  SMvertex* v2;

  if (num_e_visited == 0) // start
  {
    op_start++;

    if (edge_buffer_size)
    {
      re_conn_op->encode(rmSAJFE[last_SAJFE], SMC_START);
    }
    else
    {
      re_conn_op->encode(rmDone, 0); // more to come
    }
    last_SAJFE = SMC_START;

    PRINT_DEBUG_OUTPUT(stderr, "%d ", SMC_START);

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

    v0 = vertices[rot];
    v1 = vertices[(rot+1)%3];
    v2 = vertices[(rot+2)%3];

    // encode vertex v0

    if (v0->use_count)
    {
      // an old vertex
      re_conn->encode(rmS_0Old, 1);
      // encode its index among all active vertices
      dv_index = dv->getRelativeIndex(v0);
      re_conn_index->encode(dv->size(),dv_index);        
      start_non_manifold++;
      PRINT_DEBUG_OUTPUT(stderr, "SNM0 ");
    }
    else
    {
      // a new vertex
      re_conn->encode(rmS_0Old, 0);
      // encode its position
      compressVertexPosition(v0->v);
      // insert it into the indexable data structure
      dv->addElement(v0);
    }

    // encode vertex v1

    if (v1->use_count)
    {
      // an old vertex
      re_conn->encode(rmS_1Old, 1);
      // encode its index among all active vertices
      dv_index = dv->getRelativeIndex(v1);
      re_conn_index->encode(dv->size(),dv_index);        
      start_non_manifold++;
      PRINT_DEBUG_OUTPUT(stderr, "SNM1 ");
    }
    else
    {
      // a new vertex
      re_conn->encode(rmS_1Old, 0);
      // encode its position
      compressVertexPosition(v0->v, v1->v);
      // insert it into the indexable data structure
      dv->addElement(v1);
    }

    // encode vertex v2

    if (v2->use_count)
    {
      // an old vertex
      re_conn->encode(rmS_2Old, 1);
      // encode its index among all active vertices
      dv_index = dv->getRelativeIndex(v2);
      re_conn_index->encode(dv->size(),dv_index);        
      start_non_manifold++;
      PRINT_DEBUG_OUTPUT(stderr, "SNM2 ");
    }
    else
    {
      // a new vertex
      re_conn->encode(rmS_2Old, 0);
      // encode its position
      compressVertexPosition(v0->v, v2->v);
      // insert it into the indexable data structure
      dv->addElement(v2);
    }
  }
  else if (num_e_visited == 1) // add (or join)
  {
    op_add_join++;
    re_conn_op->encode(rmSAJFE[last_SAJFE], SMC_ADD_JOIN);
    last_SAJFE = SMC_ADD_JOIN;

    PRINT_DEBUG_OUTPUT(stderr, "%d ", SMC_ADD_JOIN);

    // rotate triangle if necessary

    rot = last_e_visited;

    // encode the (rotated) add (or join) configuration

    v0 = vertices[rot];
    v1 = vertices[(rot+1)%3];
    v2 = vertices[(rot+2)%3];

    SMedge* e0 = edges[rot];

    // encode the index of v0 (-> improve with cache of with prediction)
    lc_pos = lc->pos(v0->index);
    if (lc_pos == -1)
    {
      add_miss++;
      if (lc->pos(v1->index) != -1)
      {
        add_miss_hit_possible++;
      }
      re_conn_cache->encode(rmAJ_CacheHit,0);        
      dv_index = dv->getRelativeIndex(v0);
      re_conn_index->encode(dv->size(),dv_index);        
      PRINT_DEBUG_OUTPUT(stderr, "m%d ", dv_index);
    }
    else
    {
      add_hit++;
      add_pos[lc_pos]++;
      re_conn_cache->encode(rmAJ_CacheHit,1);        
      re_conn_cache->encode(rmAJ_CachePos,lc_pos);        
      PRINT_DEBUG_OUTPUT(stderr, "h%d ", lc_pos);
    }

    // encode which of the edges incident to v0 is e0
    candidate = -1;
    candidate_count = 0;
    if (e0->degree == 1) // is this edge only used by one other triangle?
    {
      re_conn->encode(rmAJ_Manifold,0); // is attached in a manifold way
      if (e0->target == v0) // the target vertex of this edge is v0 
      {
        re_conn->encode(rmAJ_ManifoldOriented,0); // is attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "mo ");
        // how many oriented degree-one edges are around that vertex?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->target == v0)
          {
            if (v0->list[i] == e0) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
      else
      {
        re_conn->encode(rmAJ_ManifoldOriented,1); // is attached in an not-oriented way
        PRINT_DEBUG_OUTPUT(stderr, "mno ");
        // how many not-oriented degree-one edges are around that vertex?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->target != v0)
          {
            if (v0->list[i] == e0) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
    }
    else
    {
      re_conn->encode(rmAJ_Manifold,1); // is attached in non-manifold way     
      if (e0->target == v0) // the target vertex of this edge is v0 
      {
        re_conn->encode(rmAJ_NonManifoldOriented,0); // is attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "nmo ");
        // how many oriented higher-degree edges are around that vertex?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->target == v0)
          {
            if (v0->list[i] == e0) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
      else
      {
        re_conn->encode(rmAJ_NonManifoldOriented,1); // is attached in an not-oriented way
        PRINT_DEBUG_OUTPUT(stderr, "nmno ");
        // how many not-oriented higher-degree edges are around that vertex?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->target != v0)
          {
            if (v0->list[i] == e0) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
    }
    if (candidate_count != 1)
    {
      re_conn->encode(candidate_count,candidate);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }
    // encode v2
    if (v2->use_count)
    {
      // an old vertex
      re_conn->encode(rmAJ_Old, 1);
      // encode its index among all active vertices
      dv_index = dv->getRelativeIndex(v2);
      re_conn_index->encode(dv->size(),dv_index);        
      add_non_manifold++;
      PRINT_DEBUG_OUTPUT(stderr, "AJNM %d ",dv_index);
    }
    else
    {
      // a new vertex
      re_conn->encode(rmAJ_Old, 0);
      // encode its position
      compressVertexPosition(v0->v, e0->across, v1->v, v2->v);
      // insert it into the indexable data structure
      dv->addElement(v2);
    }
  }
  else // fill or end
  {
    op_fill_end++;
     re_conn_op->encode(rmSAJFE[last_SAJFE], SMC_FILL_END);
    last_SAJFE = SMC_FILL_END;

    PRINT_DEBUG_OUTPUT(stderr, "%d ", SMC_FILL_END);

    // rotate triangle if necessary

    if (num_e_visited == 2)
    {
      rot = (last_e_unvisited + 2) % 3;
    }
    else
    {
      for (i = 0; i < 3; i++)
      {
        if (lc->pos(vertices[i]->index) != -1)
        {
          rot = i;
          break;
        }
      }
    }

    // encode the (rotated) fill (or end) configuration

    v0 = vertices[rot];
    v1 = vertices[(rot+1)%3];
    v2 = vertices[(rot+2)%3];

    SMedge* e0 = edges[rot];
//    SMedge* e1 = edges[(rot+1)%3];
    SMedge* e2 = edges[(rot+2)%3];

    // encode the index of v0
    lc_pos = lc->pos(v0->index);
    if (lc_pos == -1)
    {
      fill_miss++;
      if (lc->pos(v1->index) != -1)
      {
        fill_miss_hit_possible_v1++;
      }
      else if (lc->pos(v2->index) != -1)
      {
        fill_miss_hit_possible_v2++;
      }
      re_conn_cache->encode(rmFE_CacheHit,0);        
      dv_index = dv->getRelativeIndex(v0);
      re_conn_index->encode(dv->size(),dv_index);
      PRINT_DEBUG_OUTPUT(stderr, "m%d ", dv_index);
    }
    else
    {
      fill_hit++;
      fill_pos[lc_pos]++;
      re_conn_cache->encode(rmFE_CacheHit,1);        
      re_conn_cache->encode(rmFE_CachePos,lc_pos);        
      PRINT_DEBUG_OUTPUT(stderr, "h%d ", lc_pos);
    }

    // encode edge e0 (e.g. the edge from v0 to v1)
    candidate = -1;
    candidate_count = 0;
    if (e0->degree == 1)
    {
      re_conn->encode(rmFE_0Manifold,0); // along e0 attached in a manifold way
      if (e0->target == v0)
      {
        re_conn->encode(rmFE_0ManifoldOriented,0); // along e0 attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "mo ");
        // how many oriented manifold edges are around v0 that could potentially be e0?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->target == v0)
          {
            if (v0->list[i] == e0) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
      else
      {
        re_conn->encode(rmFE_0ManifoldOriented,1); // along e0 *not* attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "mno ");
        // how many not-oriented manifold edges are around v0 that could potentially be e0?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->target != v0)
          {
            if (v0->list[i] == e0) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
    }
    else
    {
      re_conn->encode(rmFE_0Manifold,1); // along e0 *not* attached in a manifold way
      if (e0->target == v0)
      {
        re_conn->encode(rmFE_0NonManifoldOriented,0); // along e0 attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "nmo ");
        // how many oriented non-manifold edges are around v0 that could potentially be e0?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->target == v0)
          {
            if (v0->list[i] == e0) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
      else
      {
        re_conn->encode(rmFE_0NonManifoldOriented,1); // along e0 *not* attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "nmno ");
        // how many not-oriented non-manifold edges are around v0 that could potentially be e0?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->target != v0)
          {
            if (v0->list[i] == e0) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
    }
    if (candidate_count != 1)
    {
      re_conn->encode(candidate_count,candidate);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }

    // encode edge e2 (e.g. the edge from v2 to v0)
    candidate = -1;
    candidate_count = 0;
    if (e2->degree == 1)
    {
      re_conn->encode(rmFE_2Manifold,0); // along e2 attached in a manifold way
      if (e2->origin == v0)
      {
        re_conn->encode(rmFE_2ManifoldOriented,0); // along e2 attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "mo ");
        // how many oriented manifold edges are around v0 that could potentially be e2?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->origin == v0)
          {
            if (v0->list[i] == e2) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
      else
      {
        re_conn->encode(rmFE_2ManifoldOriented,1); // along e2 *not* attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "mno ");
        // how many not-oriented manifold edges are around v0 that could potentially be e2?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->origin != v0)
          {
            if (v0->list[i] == e2) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
    }
    else
    {
      re_conn->encode(rmFE_2Manifold,1); // along e2 *not* attached in a manifold way
      if (e2->origin == v0)
      {
        re_conn->encode(rmFE_2NonManifoldOriented,0); // along e2 attached in an oriented way
        PRINT_DEBUG_OUTPUT(stderr, "nmo ");
        // how many oriented non-manifold edges are around v0 that could potentially be e2?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->origin == v0)
          {
            if (v0->list[i] == e2) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
      else
      {
        re_conn->encode(rmFE_2NonManifoldOriented,1); // along e2 *not* attached in an oriented way
        // how many not-oriented non-manifold edges are around v0 that could potentially be e2?
        PRINT_DEBUG_OUTPUT(stderr, "nmno ");
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->origin != v0)
          {
            if (v0->list[i] == e2) candidate = candidate_count;
            candidate_count++;
          }
        }
      }
    }
    if (candidate_count != 1)
    {
      re_conn->encode(candidate_count,candidate);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }
  }

  // update cache

#if defined(USE_CACHE)
  lc->put(v0->index, v1->index, v2->index, v0, v1, v2);
#endif

  if (rot != 0)
  {
//    printf("WARNING ");
  }

  // increment vertex use_counts, create edges, and update edge degrees

  for (j = 0; j < 3; j++)
  {
    i = (j+rot)%3;

    vertices[i]->use_count++;

    if (e_visited[i] == 0)
    {
      edges[i] = allocEdge(vertices[(i+2)%3]->v);
      edges[i]->origin = vertices[i];
      edges[i]->target = vertices[(i+1)%3];
      addEdgeToVertex(edges[i],vertices[i]);
      addEdgeToVertex(edges[i],vertices[(i+1)%3]);
    }
    else
    {
      if (edges[i]->degree == 1)
      {
        edges[i]->origin->degree_one--;
        edges[i]->target->degree_one--;
      }
      edges[i]->degree++;
    }
  }

  // write finalization info and finalize vertices

  for (j = 0; j < 3; j++)
  {
    i = (j+rot)%3;

    degree_one = vertices[i]->degree_one;
    if (degree_one >= MAX_DEGREE_ONE)
    {
      degree_one = MAX_DEGREE_ONE-1;
    }
    use_count = vertices[i]->use_count;
    if (use_count >= MAX_USE_COUNT)
    {
      use_count = MAX_USE_COUNT-1;
    }
    re_conn_final->encode(rmFinalized[degree_one][use_count],t_final[i]);

//    PRINT_DEBUG_OUTPUT(stderr, "f%d:%d:%d ",degree_one,use_count,t_final[i]);

    if (t_final[i])
    {
      SMvertex* vertex = vertices[i];
      vertex_hash->erase(hash_elements[i]);
      dv->removeElement(vertex);
      for (int k = 0; k < vertex->list_size; k++)
      {
        SMedge* edge = vertex->list[k];
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

void SMwriter_smc_old::write_finalized(int final_idx)
{
  fprintf(stderr, "ERROR: write_finalized(int final_idx) not supported by SMwriter_smc_old\n");
  exit(0);
}

bool SMwriter_smc_old::open(FILE* file, int nbits)
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

  initVertexBuffer(1024);
  initEdgeBuffer(1024);
  vertex_hash = new my_vertex_hash;
  dv = new DynamicVector();
  lc = new LittleCache();

  initEncoder(file);
  initModels(1);

  // write precision
  re_conn->encode(25,nbits);

  // store precision (is needed in write_header)
  precision_bits = nbits;

  v_count = 0;
  f_count = 0;

  return true;
}

void SMwriter_smc_old::close(bool close_file, bool update_header)
{
  // close of SMwriter_smc_old
  re_conn_op->encode(rmDone, 1); // done
  PRINT_CONTROL_OUTPUT(stderr,"edge_buffer_size %d edge_buffer_maxsize %d\n",edge_buffer_size,edge_buffer_maxsize);
  PRINT_CONTROL_OUTPUT(stderr,"vertex_buffer_size %d vertex_buffer_maxsize %d\n",vertex_buffer_size,vertex_buffer_maxsize);
  PRINT_CONTROL_OUTPUT(stderr,"op_start %d op_add_join %d op_fill_end %d\n",op_start,op_add_join,op_fill_end);
  PRINT_CONTROL_OUTPUT(stderr,"start_non_manifold %d add_non_manifold %d\n",start_non_manifold,add_non_manifold);
  PRINT_CONTROL_OUTPUT(stderr,"add_miss %d (hit_possible %d) add_hit %d (%d %d %d)\n",add_miss,add_miss_hit_possible,add_hit,add_pos[0],add_pos[1],add_pos[2]);
  PRINT_CONTROL_OUTPUT(stderr,"fill_miss %d (hit_possible %d %d) fill_hit %d  (%d %d %d)\n",fill_miss,fill_miss_hit_possible_v1,fill_miss_hit_possible_v2,fill_hit,fill_pos[0],fill_pos[1],fill_pos[2]);
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

  if (dv->size()) fprintf(stderr,"WARNING: there are %d unfinalized vertices\n         this mesh will currently not decompress correctly\n",dv->size());
  if (vertex_hash->size()-dv->size()) fprintf(stderr,"WARNING: there are %d unused vertices\n         these vertices have not been compressed\n",vertex_hash->size()-dv->size());

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

void SMwriter_smc_old::write_header()
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
    re_conn->encode(2,0);
  }
  else
  {
    re_conn->encode(2,1);
    re_conn->encodeFloat(bb_min_f[0]);
    re_conn->encodeFloat(bb_min_f[1]);
    re_conn->encodeFloat(bb_min_f[2]);
  }
  if (bb_max_f == 0)
  {
    re_conn->encode(2,0);
  }
  else
  {
    re_conn->encode(2,1);
    re_conn->encodeFloat(bb_max_f[0]);
    re_conn->encodeFloat(bb_max_f[1]);
    re_conn->encodeFloat(bb_max_f[2]);
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

  if (bb_min_f && bb_max_f && precision_bits)
  {
    if (true) // if want to use integer quantization
    {
      re_conn->encode(2,0);
    }
    else
    {
      /* yet to be implemented */
    }

    pq = new PositionQuantizer();
    pq->SetMinMax(bb_min_f, bb_max_f);
    pq->SetPrecision(precision_bits);
    pq->SetupQuantizer();

    ic[0] = new IntegerCompressor();
    ic[1] = new IntegerCompressor();
    ic[2] = new IntegerCompressor();

    ic[0]->SetMax(pq->m_aiRangeCode[0], pq->m_aiAbsRangeCorrector[0], pq->m_aiAbsRangeCorrector[0]);
    ic[1]->SetMax(pq->m_aiRangeCode[1], pq->m_aiAbsRangeCorrector[1], pq->m_aiAbsRangeCorrector[1]);
    ic[2]->SetMax(pq->m_aiRangeCode[2], pq->m_aiAbsRangeCorrector[2], pq->m_aiAbsRangeCorrector[2]);

    ic[0]->SetPrecision(precision_bits);
    ic[0]->SetupCompressor(re_geom);
    ic[1]->SetPrecision(precision_bits);
    ic[2]->SetPrecision(precision_bits);
    ic[1]->SetupCompressor(re_geom);
    ic[2]->SetupCompressor(re_geom);
  }
  else
  {
    pq = 0;
    fc[0] = new FloatCompressor();
    fc[1] = new FloatCompressor();
    fc[2] = new FloatCompressor();

    fc[0]->SetPrecision(precision_bits);
    fc[0]->SetupCompressor(re_geom,0);
    fc[1]->SetPrecision(precision_bits);
    fc[2]->SetPrecision(precision_bits);
    fc[1]->SetupCompressor(re_geom,0);
    fc[2]->SetupCompressor(re_geom,0);
  }
}

SMwriter_smc_old::SMwriter_smc_old()
{
  // re-init of SMwriter interface
  need_pre_order = true;
  // init of SMwriter_smc_old interface
}

SMwriter_smc_old::~SMwriter_smc_old()
{
  // clean-up for SMwriter_smc_old
}
