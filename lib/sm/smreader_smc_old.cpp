/*
===============================================================================

  FILE:  smreader_smc_old.cpp
  
  CONTENTS:
  
    Reads a mesh from our compressed Streaming Mesh format (SMC) and provides
    access to it in form of a Streaming Mesh.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    07 October 2003 -- initial version created the day of the California recall
  
===============================================================================
*/
#include "smreader_smc_old.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "vec3fv.h"
#include "vec3iv.h"

#include "rangemodel.h"
#include "rangedecoder.h"

#include "dynamicvector.h"
#include "littlecache.h"

#include "floatcompressor.h"

#include "positionquantizer.h"
#include "integercompressor.h"

//#define PRINT_CONTROL_OUTPUT fprintf
#define PRINT_CONTROL_OUTPUT if (0) fprintf
//#define PRINT_DEBUG_OUTPUT fprintf
#define PRINT_DEBUG_OUTPUT if (0) fprintf

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
  SMedge* buffer_next; // used for efficient memory management
  SMvertex* origin;
  SMvertex* target;
  int degree;
  float across[3];
} SMedge;

static DynamicVector* dv;

static LittleCache* lc;

static FloatCompressor* fc[3];

static PositionQuantizer* pq;
static IntegerCompressor* ic[3];

// rangecoder and probability tables

#define MAX_DEGREE_ONE 4
#define MAX_USE_COUNT 15

static RangeDecoder* rd_conn;
static RangeDecoder* rd_conn_op;
static RangeDecoder* rd_conn_cache;
static RangeDecoder* rd_conn_index;
static RangeDecoder* rd_conn_final;

static RangeDecoder* rd_geom;

// is there more to encode
static RangeModel* rmDone;

// which operation
static int last_SAJFE = 0;
static RangeModel** rmSAJFE;

// used for start operations
static RangeModel* rmS_0Old;
static RangeModel* rmS_1Old;
static RangeModel* rmS_2Old;

//static RangeModel* rmS_CacheHit;
//static RangeModel* rmS_CachePos;

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

static void initDecoder(FILE* file)
{
  if (file == 0)
  {
    fprintf(stderr,"FATAL ERROR: file pointer is zero\n");
    exit(0);
  }
  rd_conn = new RangeDecoder(file);
  rd_conn_op = rd_conn;
  rd_conn_cache = rd_conn;
  rd_conn_index = rd_conn;
  rd_conn_final = rd_conn;
  rd_geom = rd_conn;
}

static void finishDecoder()
{
  rd_conn->done();
  delete rd_conn;
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
  int i,j;

  // range table for whether there is more
  delete rmDone;

  // range table for which operation it is
  delete rmSAJFE[0];
  delete rmSAJFE[1];
  delete rmSAJFE[2];
  free(rmSAJFE);

  // range tables for start (S) operation
  delete rmS_0Old;
  delete rmS_1Old;
  delete rmS_2Old;

  // range tables for add (A) / join (J) operations
  delete rmAJ_Manifold;
  delete rmAJ_ManifoldOriented;
  delete rmAJ_NonManifoldOriented;
  delete rmAJ_Old;

  delete rmAJ_CacheHit;
  delete rmAJ_CachePos;

  // range tables for fill (F) / end (E) operations
  delete rmFE_0Manifold;
  delete rmFE_0ManifoldOriented;
  delete rmFE_0NonManifoldOriented;
  delete rmFE_2Manifold;
  delete rmFE_2ManifoldOriented;
  delete rmFE_2NonManifoldOriented;

  delete rmFE_CacheHit;
  delete rmFE_CachePos;

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

static int op_start = 0;
static int op_add_join = 0;
static int op_fill_end = 0;

static int add_hit = 0;
static int add_miss = 0;
static int add_pos[] = {0,0,0};

static int fill_hit = 0;
static int fill_miss = 0;
static int fill_pos[] = {0,0,0};

static int start_non_manifold = 0;
static int add_non_manifold = 0;

static int vertex_buffer_size = 0;
static int edge_buffer_size = 0;

static int vertex_buffer_maxsize = 0;
static int edge_buffer_maxsize = 0;

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

static void decompressVertexPosition(float* n)
{
  PRINT_CONTROL_OUTPUT(stderr, "none\n");
  if (pq)
  {
    int* qn = (int*)n;
    for (int i = 0; i < 3; i++)
    {
      qn[i] = ic[i]->DecompressNone();
    }
  }
  else
  {
    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->DecompressNone();
    }
    PRINT_DEBUG_OUTPUT(stderr, "n %f %f %f\n",n[0],n[1],n[2]);
  }
}

static void decompressVertexPosition(float* l, float* n)
{
  PRINT_DEBUG_OUTPUT(stderr, "last\n");
  PRINT_DEBUG_OUTPUT(stderr, "l %f %f %f\n",l[0],l[1],l[2]);
  if (pq)
  {
    int* ql = (int*)l;
    int* qn = (int*)n;
    int corr[3];
    for (int i = 0; i < 3; i++)
    {
      corr[i] = ic[i]->DecompressLast();
    }
    pq->AddCorrector(ql, corr, qn);
  }
  else
  {
    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->DecompressLast(l[i]);
    }
    PRINT_DEBUG_OUTPUT(stderr, "n %f %f %f\n",n[0],n[1],n[2]);
  }
}

static void decompressVertexPosition(float* a, float* b, float* c, float* n)
{
  if (pq)
  {
    int* qa = (int*)a;
    int* qb = (int*)b;
    int* qc = (int*)c;
    int* qn = (int*)n;
    int pred[3];
    VecAdd3iv(pred, qa, qc);
    VecSelfSubtract3iv(pred, qb);
    int corr[3];
    for (int i = 0; i < 3; i++)
    {
      corr[i] = ic[i]->DecompressAcross();
    }
    pq->AddCorrector(pred, corr, qn); 
  }
  else
  {
    float pp[3];
    VecAdd3fv(pp, a, c);
    VecSelfSubtract3fv(pp, b);
    PRINT_DEBUG_OUTPUT(stderr, "across\n");
    PRINT_DEBUG_OUTPUT(stderr, "a %f %f %f\n",a[0],a[1],a[2]);
    PRINT_DEBUG_OUTPUT(stderr, "b %f %f %f\n",b[0],b[1],b[2]);
    PRINT_DEBUG_OUTPUT(stderr, "c %f %f %f\n",c[0],c[1],c[2]);
    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->DecompressAcross(pp[i]);
    }
    PRINT_DEBUG_OUTPUT(stderr, "n %f %f %f\n",n[0],n[1],n[2]);
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

bool SMreader_smc_old::open(FILE* file)
{
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

  initVertexBuffer(1024);
  initEdgeBuffer(1024);
  dv = new DynamicVector();
  lc = new LittleCache();
  initDecoder(file);
  initModels(0);

  // read precision
  nbits = rd_conn->decode(25);

  v_count = 0;
  f_count = 0;

  read_header();

  return true;
}

void SMreader_smc_old::close(bool close_file)
{
  // close of SMreader_smc_old
  finishDecoder();
  finishModels();
  if (pq)
  {
    ic[0]->FinishDecompressor();
    ic[1]->FinishDecompressor();
    ic[2]->FinishDecompressor();
    delete ic[0];
    delete ic[1];
    delete ic[2];
    delete pq;
  }
  else
  {
    delete fc[0];
    delete fc[1];
    delete fc[2];
  }

  if (dv->size()) fprintf(stderr,"WARNING: there are %d unfinalized vertices\n",dv->size());

  delete dv;
  delete lc;

  if (close_file && file /*&& file != stdin*/) fclose(file);

  nbits = -1;
  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;

  PRINT_CONTROL_OUTPUT(stderr,"edge_buffer_size %d edge_buffer_maxsize %d\n",edge_buffer_size,edge_buffer_maxsize);
  PRINT_CONTROL_OUTPUT(stderr,"vertex_buffer_size %d vertex_buffer_maxsize %d\n",vertex_buffer_size,vertex_buffer_maxsize);
  PRINT_CONTROL_OUTPUT(stderr,"op_start %d op_add_join %d op_fill_end %d\n",op_start,op_add_join,op_fill_end);
  PRINT_CONTROL_OUTPUT(stderr,"start_non_manifold %d add_non_manifold %d\n",start_non_manifold,add_non_manifold);
  PRINT_CONTROL_OUTPUT(stderr,"add_miss %d (hit_possible --) add_hit %d (%d %d %d)\n",add_miss,add_hit,add_pos[0],add_pos[1],add_pos[2]);
  PRINT_CONTROL_OUTPUT(stderr,"fill_miss %d (hit_possible -- --) fill_hit %d  (%d %d %d)\n",fill_miss,fill_hit,fill_pos[0],fill_pos[1],fill_pos[2]);

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

void SMreader_smc_old::read_header()
{
  // read nverts
  if (rd_conn->decode(2))
  {
    nverts = rd_conn->decodeInt();
    PRINT_CONTROL_OUTPUT(stderr,"nverts: %d\n",nverts);
  }

  // read nfaces
  if (rd_conn->decode(2))
  {
    nfaces = rd_conn->decodeInt();
    PRINT_CONTROL_OUTPUT(stderr,"nfaces: %d\n",nfaces);
  }

  bool has_bb = true;

  // read bounding box
  if (rd_conn->decode(2))
  {
    if (bb_min_f == 0) bb_min_f = new float[3];
    bb_min_f[0] = rd_conn->decodeFloat();
    bb_min_f[1] = rd_conn->decodeFloat();
    bb_min_f[2] = rd_conn->decodeFloat();
  }
  else
  {
    // no bounding box min
    has_bb = false;
  }
  if (rd_conn->decode(2))
  {
    if (bb_max_f == 0) bb_max_f = new float[3];
    bb_max_f[0] = rd_conn->decodeFloat();
    bb_max_f[1] = rd_conn->decodeFloat();
    bb_max_f[2] = rd_conn->decodeFloat();
  }
  else
  {
    // no bounding box max
    has_bb = false;
  }

  // read comments
  if (rd_conn->decode(2))
  {
    /* yet to be implemented */
  }

  // setup geometry compressor

  if (has_bb && nbits)
  {
    if (rd_conn->decode(2) == 0) // if want to use integer quantization
    {
      pq = new PositionQuantizer();

      pq->SetMinMax(bb_min_f, bb_max_f);
      pq->SetPrecision(nbits);
      pq->SetupQuantizer();

      ic[0] = new IntegerCompressor();
      ic[1] = new IntegerCompressor();
      ic[2] = new IntegerCompressor();

      ic[0]->SetMax(pq->m_aiRangeCode[0], pq->m_aiAbsRangeCorrector[0], pq->m_aiAbsRangeCorrector[0]);
      ic[1]->SetMax(pq->m_aiRangeCode[1], pq->m_aiAbsRangeCorrector[1], pq->m_aiAbsRangeCorrector[1]);
      ic[2]->SetMax(pq->m_aiRangeCode[2], pq->m_aiAbsRangeCorrector[2], pq->m_aiAbsRangeCorrector[2]);

      ic[0]->SetPrecision(nbits);
      ic[1]->SetPrecision(nbits);
      ic[2]->SetPrecision(nbits);

      ic[0]->SetupDecompressor(rd_geom);
      ic[1]->SetupDecompressor(rd_geom);
      ic[2]->SetupDecompressor(rd_geom);
    }
    else
    {
      /* yet to be implemented */
    }
  }
  else
  {
    pq = 0;
    fc[0] = new FloatCompressor();
    fc[1] = new FloatCompressor();
    fc[2] = new FloatCompressor();

    fc[0]->SetPrecision(nbits);
    fc[1]->SetPrecision(nbits);
    fc[2]->SetPrecision(nbits);

    fc[0]->SetupDecompressor(rd_geom,0);
    fc[1]->SetupDecompressor(rd_geom,0);
    fc[2]->SetupDecompressor(rd_geom,0);
  }
}

int SMreader_smc_old::decompress_triangle()
{
  int i,j;
  int op;
  int lc_pos;
  int dv_index;
  int degree_one;
  int use_count;

  int candidate;
  int candidate_count;
  SMedge* candidates[100];

  SMvertex* vertices[3];
  SMedge* edges[3];

  SMvertex* v0;
  SMvertex* v1;
  SMvertex* v2;

  SMedge* e0;
//  SMedge* e1;
  SMedge* e2;

  if (edge_buffer_size)
  {
    op = rd_conn_op->decode(rmSAJFE[last_SAJFE]);
  }
  else
  {
    PRINT_DEBUG_OUTPUT(stderr, "ZERO ");
    if (rd_conn_op->decode(rmDone))
    {
      return 0;
    }
    else
    {
      PRINT_DEBUG_OUTPUT(stderr, "GO ");
      op = SMC_START;
    }
  }
  last_SAJFE = op;

  have_triangle = 1;
  
  
  PRINT_DEBUG_OUTPUT(stderr, "%d ", op);

  if (op == SMC_ADD_JOIN)
  {
    op_add_join++;

    // encode the index of v0 (-> improve with cache of with prediction)
    if (rd_conn_cache->decode(rmAJ_CacheHit) == 0)
    {
      add_miss++;
      dv_index = rd_conn_index->decode(dv->size());        
      v0 = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      PRINT_DEBUG_OUTPUT(stderr, "m%d ", dv_index);
    }
    else
    {
      add_hit++;
      lc_pos = rd_conn_cache->decode(rmAJ_CachePos);
      add_pos[lc_pos]++;
      v0 = (SMvertex*)lc->get(lc_pos);
      PRINT_DEBUG_OUTPUT(stderr, "h%d ", lc_pos);
    }

    // decode which of the edges incident to v0 is e0
    candidate = 0;
    candidate_count = 0;
    if (rd_conn->decode(rmAJ_Manifold) == 0) // is this edge only used by one other triangle?
    {
      if (rd_conn->decode(rmAJ_ManifoldOriented) == 0) // is attached in an oriented way
      {
        PRINT_DEBUG_OUTPUT(stderr, "mo ");
        // how many oriented degree-one edges are around that vertex?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->target == v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
      else
      {
        PRINT_DEBUG_OUTPUT(stderr, "mno ");
        // how many not-oriented degree-one edges are around that vertex?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->target != v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
    }
    else
    {
      if (rd_conn->decode(rmAJ_NonManifoldOriented) == 0) // is attached in an oriented way
      {
        PRINT_DEBUG_OUTPUT(stderr, "nmo ");
        // how many oriented higher-degree edges are around that vertex?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->target == v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
      else
      {
        PRINT_DEBUG_OUTPUT(stderr, "nmno ");
        // how many not-oriented higher-degree edges are around that vertex?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->target != v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
    }
    if (candidate_count != 1)
    {
      candidate = rd_conn->decode(candidate_count);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }
    e0 = candidates[candidate];
    v1 = (e0->target == v0 ? e0->origin : e0->target);

    // encode v2
    if (rd_conn->decode(rmAJ_Old))
    {
      // an old vertex
      dv_index = rd_conn_index->decode(dv->size());        
      v2 = (SMvertex*) dv->getElementWithRelativeIndex(dv_index);
      add_non_manifold++;
      PRINT_DEBUG_OUTPUT(stderr, "AJNM %d ",dv_index);
    }
    else
    {
      // a new vertex
      v2 = allocVertex();
      // give it its index
      v2->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(v0->v, e0->across, v1->v, v2->v);
      // insert it into the indexable data structure
      dv->addElement(v2);
      new_vertices[have_new] = 2;
      have_new++;
    }
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    edges[0] = e0;
    edges[1] = 0;
    edges[2] = 0;
  }
  else if (op == SMC_FILL_END)
  {
    op_fill_end++;

    // encode the index of v0
    if (rd_conn_cache->decode(rmFE_CacheHit) == 0)
    {
      fill_miss++;
      dv_index = rd_conn_index->decode(dv->size());
      v0 = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      PRINT_DEBUG_OUTPUT(stderr, "m%d ", dv_index);
    }
    else
    {
      fill_hit++;
      lc_pos = rd_conn_cache->decode(rmFE_CachePos);        
      fill_pos[lc_pos]++;
      v0 = (SMvertex*)lc->get(lc_pos);
      PRINT_DEBUG_OUTPUT(stderr, "h%d ", lc_pos);
    }

    // decode edge e0 (e.g. the edge from v0 to v1)
    candidate = 0;
    candidate_count = 0;
    if (rd_conn->decode(rmFE_0Manifold) == 0)  // along e0 attached in a manifold way
    {
      if (rd_conn->decode(rmFE_0ManifoldOriented) == 0) // along e0 attached in an oriented way
      {
        PRINT_DEBUG_OUTPUT(stderr, "mo ");
        // how many oriented manifold edges are around v0 that could potentially be e0?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->target == v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
      else
      {
        PRINT_DEBUG_OUTPUT(stderr, "mno ");
        // how many not-oriented manifold edges are around v0 that could potentially be e0?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->target != v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
    }
    else
    {
      if (rd_conn->decode(rmFE_0NonManifoldOriented) == 0) // along e0 attached in an oriented way
      {
        PRINT_DEBUG_OUTPUT(stderr, "nmo ");
        // how many oriented non-manifold edges are around v0 that could potentially be e0?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->target == v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
      else
      {
        PRINT_DEBUG_OUTPUT(stderr, "nmno ");
        // how many not-oriented non-manifold edges are around v0 that could potentially be e0?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->target != v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
    }
    if (candidate_count != 1)
    {
      candidate = rd_conn->decode(candidate_count);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }
    e0 = candidates[candidate];
    v1 = (e0->target == v0 ? e0->origin : e0->target);

    // encode edge e2 (e.g. the edge from v2 to v0)
    candidate = 0;
    candidate_count = 0;
    if (rd_conn->decode(rmFE_2Manifold) == 0) // along e2 attached in a manifold way
    {
      if (rd_conn->decode(rmFE_2ManifoldOriented) == 0) // along e2 attached in an oriented way
      {
        PRINT_DEBUG_OUTPUT(stderr, "mo ");
        // how many oriented manifold edges are around v0 that could potentially be e2?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->origin == v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
      else
      {
        PRINT_DEBUG_OUTPUT(stderr, "mno ");
        // how many not-oriented manifold edges are around v0 that could potentially be e2?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree == 1 && v0->list[i]->origin != v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
    }
    else
    {
      if (rd_conn->decode(rmFE_2NonManifoldOriented) == 0) // along e2 attached in an oriented way
      {
        PRINT_DEBUG_OUTPUT(stderr, "nmo ");
        // how many oriented non-manifold edges are around v0 that could potentially be e2?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->origin == v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
      else
      {
        PRINT_DEBUG_OUTPUT(stderr, "nmno ");
        // how many not-oriented non-manifold edges are around v0 that could potentially be e2?
        for (i = 0; i < v0->list_size; i++)
        {
          if (v0->list[i]->degree > 1 && v0->list[i]->origin != v0)
          {
            candidates[candidate_count] = v0->list[i];
            candidate_count++;
          }
        }
      }
    }
    if (candidate_count != 1)
    {
      candidate = rd_conn->decode(candidate_count);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
    }
    e2 = candidates[candidate];
    v2 = (e2->origin == v0 ? e2->target : e2->origin);

    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    edges[0] = e0;
    edges[1] = 0;
    edges[2] = e2;

    // is it a FILL or an END operation

    for (i = 0; i < v2->list_size; i++)
    {
       if ((v2->list[i]->target == v1) || (v2->list[i]->origin == v1))
       {
         edges[1] = v2->list[i];
       }
    }
  }
  else if (op == SMC_START)
  {
    op_start++;

    // decode vertex v0

    if (rd_conn->decode(rmS_0Old))
    {
      // an old vertex
      dv_index = rd_conn_index->decode(dv->size());        
      v0 = (SMvertex*) dv->getElementWithRelativeIndex(dv_index);
      start_non_manifold++;
      PRINT_DEBUG_OUTPUT(stderr, "SNM0 ");
    }
    else
    {
      // a new vertex
      v0 = allocVertex();
      // give it its index
      v0->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(v0->v);
      // insert it into the indexable data structure
      dv->addElement(v0);
      new_vertices[have_new] = 0;
      have_new++;
    }

    // decode vertex v1

    if (rd_conn->decode(rmS_1Old))
    {
      // an old vertex
      dv_index = rd_conn_index->decode(dv->size());        
      v1 = (SMvertex*) dv->getElementWithRelativeIndex(dv_index);
      start_non_manifold++;
      PRINT_DEBUG_OUTPUT(stderr, "SNM1 ");
    }
    else
    {
      // a new vertex
      v1 = allocVertex();
      // give it its index
      v1->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(v0->v, v1->v);
      // insert it into the indexable data structure
      dv->addElement(v1);
      new_vertices[have_new] = 1;
      have_new++;
    }

    // decode vertex v2

    if (rd_conn->decode(rmS_2Old))
    {
      // an old vertex
      dv_index = rd_conn_index->decode(dv->size());        
      v2 = (SMvertex*) dv->getElementWithRelativeIndex(dv_index);
      start_non_manifold++;
      PRINT_DEBUG_OUTPUT(stderr, "SNM2 ");
    }
    else
    {
      // a new vertex
      v2 = allocVertex();
      // give it its index
      v2->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(v0->v, v2->v);
      // insert it into the indexable data structure
      dv->addElement(v2);
      new_vertices[have_new] = 2;
      have_new++;
    }
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    edges[0] = 0;
    edges[1] = 0;
    edges[2] = 0;
  }

  #if defined(USE_CACHE)
  lc->put(v0, v1, v2);
  #endif

  // copy vertex and triangle data into API

  for (i = 0; i < 3; i++)
  {
    t_idx[i] = vertices[i]->index;
    t_pos_f[i] = vertices[i]->v;
  }

  // increment vertex use_counts, create edges, and update edge degrees

  for (i = 0; i < 3; i++)
  {
    vertices[i]->use_count++;

    if (edges[i] == 0)
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

  // read finalization info and finalize vertices

  for (i = 0; i < 3; i++)
  {
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
    
    if (rd_conn_final->decode(rmFinalized[degree_one][use_count]))
    {
      t_final[i] = true;
      SMvertex* vertex = vertices[i];
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
      finalized_vertices[have_finalized] = i;
      have_finalized++;
    }
    else
    {
      t_final[i] = false;
    }
    PRINT_DEBUG_OUTPUT(stderr, "f%d:%d:%d ",degree_one,use_count,t_final[i]);
  }

  return 1;
}

SMevent SMreader_smc_old::read_element()
{
  if ((have_new + have_triangle) == 0)
  {
    next_new = 0;
    have_finalized = next_finalized = 0;
    if (decompress_triangle() == 0)
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
      return SM_EOF;
    }
  }
  if (have_new)
  {
    v_idx = t_idx[new_vertices[next_new]];
    if (pq)
    {
      pq->DeQuantize((int*) t_pos_f[new_vertices[next_new]], v_pos_f);
    }
    else
    {
      VecCopy3fv(v_pos_f, t_pos_f[new_vertices[next_new]]);
    }
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
  return SM_ERROR;
}

SMevent SMreader_smc_old::read_event()
{
  if ((have_new + have_triangle + have_finalized) == 0)
  {
    next_new = 0;
    next_finalized = 0;
    if (decompress_triangle() == 0)
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
      return SM_EOF;
    }
  }
  if (have_new)
  {
    v_idx = t_idx[new_vertices[next_new]];
    if (pq)
    {
      pq->DeQuantize((int*) t_pos_f[new_vertices[next_new]], v_pos_f);
    }
    else
    {
      VecCopy3fv(v_pos_f, t_pos_f[new_vertices[next_new]]);
    }
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

SMreader_smc_old::SMreader_smc_old()
{
  // init of SMreader_smc_old
  nbits = -1;
  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;
}

SMreader_smc_old::~SMreader_smc_old()
{
  // clean-up for SMreader_smc_old
}
