/*
===============================================================================

  FILE:  smreader_smd.cpp
  
  CONTENTS:
  
    see header file

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see header file
  
===============================================================================
*/
#include "smreader_smd.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "dynamicvector.h"
#include "dynamicqueue.h"
#include "littlecache.h"

#include "rangemodel.h"
#include "rangedecoder.h"

#include "floatcompressor.h"
#include "integercompressor_new.h"
#include "positionquantizer_new.h"

#include "vec3fv.h"
#include "vec3iv.h"

#define PRINT_CONTROL_OUTPUT
#undef PRINT_CONTROL_OUTPUT
#define PRINT_DEBUG_OUTPUT if (0) fprintf

#define USE_CACHE 1

#define SMC_ADD 0
#define SMC_JOIN 1
#define SMC_FILL_END 2
#define SMC_SKIP 3
#define SMC_BORDER 4
#define SMC_START 3

// data structures used for streaming compression

struct SMedge;

typedef struct SMvertex
{
  SMvertex* buffer_next;    // used for efficient memory management (and also for the dynamic vector)
  float v[3];
  int index;
  int use_count;
  int list_size;
  int list_alloc;
  SMedge** list;
} SMvertex;

typedef struct SMedge
{
  SMedge* buffer_next; // used for efficient memory management (and also for the dynamic queue)
  SMvertex* origin;
  SMvertex* target;
  float across[3];
} SMedge;

static int next_waiting = 100;
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

static RangeDecoder* rd_conn;
static RangeDecoder* rd_conn_op;
static RangeDecoder* rd_conn_rl;
static RangeDecoder* rd_conn_index;
static RangeDecoder* rd_conn_final;

static RangeDecoder* rd_geom;

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

static void initDecoder(FILE* file)
{
  if (file == 0)
  {
    fprintf(stderr,"FATAL ERROR: file pointer is zero\n");
    exit(0);
  }
  rd_conn = new RangeDecoder(file);
  rd_conn_op = rd_conn;
  rd_conn_rl = rd_conn;
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

  // range tables for fill (F) / end (E) operations
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

  // range table for whether there is more
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

  // range tables for start (S) operation
  delete rmWhichStart;

  // range tables for fill (F) / end (E) operations
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

static int prediction_none = 0;
static int prediction_last = 0;
static int prediction_across = 0;

static int right_confirm = 0;
static int right_correct = 0;
static int left_confirm = 0;
static int left_correct = 0;

static int vertex_buffer_maxsize = 0;
static int edge_buffer_maxsize = 0;

#endif

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

static void decompressVertexPosition(float* n)
{
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
  }
}

static void decompressVertexPosition(const float* l, float* n)
{
  if (pq)
  {
    const int* ql = (const int*)l;
    int* qn = (int*)n;
    for (int i = 0; i < 3; i++)
    {
      qn[i] = ic[i]->DecompressLast(ql[i]);
    }
  }
  else
  {
    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->DecompressLast(l[i]);
    }
  }
}

static void decompressVertexPosition(const float* a, const float* b, const float* c, float* n)
{
  if (pq)
  {
    const int* qa = (const int*)a;
    const int* qb = (const int*)b;
    const int* qc = (const int*)c;
    int* qn = (int*)n;
    int pred[3];
    VecAdd3iv(pred, qa, qc);
    VecSelfSubtract3iv(pred, qb);
    pq->Clamp(pred);
    for (int i = 0; i < 3; i++)
    {
      qn[i] = ic[i]->DecompressAcross(pred[i]);
    }
  }
  else
  {
    float pred[3];
    VecParallelogram3fv(pred, a, b, c);

    /*
    This is CRAZY Micosoft crap.... the floating-point calculation below can have different results in Debug and Release mode!!!!
    VecAdd3fv(pp, a, c);
    VecSelfSubtract3fv(pred, b);
    */

    for (int i = 0; i < 3; i++)
    {
      n[i] = fc[i]->DecompressAcross(pred[i]);
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

#define SM_VERSION 2 // this is SMD

bool SMreader_smd::open(FILE* file)
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

  int input = fgetc(file);
  // read version
  if (input != SM_VERSION)
  {
    fprintf(stderr,"ERROR: wrong SMreader (need %d but this is SMreader_smd %d)\n",input,SM_VERSION);
    if (input == EOF) fprintf(stderr,"       check ... maybe the file (or the stream) is empty\n");
    exit(0);
  }

  initVertexBuffer(1024);
  initEdgeBuffer(1024);

  traversal_queue = new DynamicQueue();
  little_cache = new LittleCache();

  dv = new DynamicVector();

  initDecoder(file);
  initModels(0);

  // read precision
  nbits = rd_conn->decode(25);

  v_count = 0;
  f_count = 0;

  read_header();

  return true;
}

void SMreader_smd::close(bool close_file)
{
  // close of SMreader_smd
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

  if (dv->size()) fprintf(stderr,"WARNING: there are %d unfinalized vertices\n         this mesh will not decompress correctly",dv->size());

  delete dv;
  delete traversal_queue;
  delete little_cache;

  if (close_file && file /*&& file != stdin*/) fclose(file);

  nbits = -1;
  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"edge_buffer_size %d edge_buffer_maxsize %d\n",edge_buffer_size,edge_buffer_maxsize);
  fprintf(stderr,"vertex_buffer_size %d vertex_buffer_maxsize %d\n",vertex_buffer_size,vertex_buffer_maxsize);
  fprintf(stderr,"start %d add %d join %d fill %d end %d skip %d border %d\n",op_start,op_add,op_join,op_fill,op_end,op_skip,op_border);
#endif

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

void SMreader_smd::read_header()
{
  // read nverts
  if (rd_conn->decode(2))
  {
    nverts = rd_conn->decodeInt();
#ifdef PRINT_CONTROL_OUTPUT
    fprintf(stderr,"nverts: %d\n",nverts);
#endif
  }

  // read nfaces
  if (rd_conn->decode(2))
  {
    nfaces = rd_conn->decodeInt();
#ifdef PRINT_CONTROL_OUTPUT
    fprintf(stderr,"nfaces: %d\n",nfaces);
#endif
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
      pq = new PositionQuantizerNew();

      pq->SetMinMax(bb_min_f, bb_max_f);
      pq->SetPrecision(nbits);
      pq->SetupQuantizer();

      ic[0] = new IntegerCompressorNew();
      ic[1] = new IntegerCompressorNew();
      ic[2] = new IntegerCompressorNew();

      ic[0]->SetRange(pq->m_aiQuantRange[0]);
      ic[1]->SetRange(pq->m_aiQuantRange[1]);
      ic[2]->SetRange(pq->m_aiQuantRange[2]);

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

bool SMreader_smd::decompress_triangle_waiting()
{
  int i, j;
  int op;
  int dv_index;
  int candidate;
  int candidate_count;

  SMedge* edges[3];
  SMvertex* vertices[3];

  if (dv->size() == 0)
  {
    if (rd_conn_op->decode(rmDone))
    {
      return false;
    }
    op = SMC_START;
  }
  else
  {
    op = rd_conn_op->decode(rmWaitingOp);
  }

  PRINT_DEBUG_OUTPUT(stderr, "w%d ", SMC_START);

  if (op == SMC_START)
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_start++;
#endif

    candidate = rd_conn_op->decode(rmWhichStart);

    // decode vertex v0

    if (candidate)
    {
      // decode v0's relative index 
      dv_index = rd_conn_index->decode(dv->size());
      // get v0 among all active vertices
      vertices[0] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      // decrement non-manifold start vertex counter
      candidate--;
    }
    else
    {
      // allocate new vertex
      vertices[0] = allocVertex();
      // give it its index
      vertices[0]->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(vertices[0]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[0]);
      new_vertices[have_new] = 0;
      have_new++;
    }

    // decode vertex v1

    if (candidate)
    {
      // decode v1's relative index 
      dv_index = rd_conn_index->decode(dv->size());
      // get v1 among all active vertices
      vertices[1] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      // decrement non-manifold start vertex counter
      candidate--;
    }
    else
    {
      // allocate new vertex
      vertices[1] = allocVertex();
      // give it its index
      vertices[1]->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(vertices[0]->v, vertices[1]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[1]);
      new_vertices[have_new] = 1;
      have_new++;
    }

    // decode vertex v2

    if (candidate)
    {
      // decode v2's relative index 
      dv_index = rd_conn_index->decode(dv->size());
      // get v1 among all active vertices
      vertices[2] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      // statistics
      candidate--;
    }
    else
    {
      // allocate new vertex
      vertices[2] = allocVertex();
      // give it its index
      vertices[2]->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(vertices[0]->v, vertices[2]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[2]);
      new_vertices[have_new] = 2;
      have_new++;
    }

    // no edge exists yet

    edges[0] = 0;
    edges[1] = 0;
    edges[2] = 0;
  }
  else
  {
    // decode v0

    // decode v0's relative index 
    dv_index = rd_conn_index->decode(dv->size());
    // get v0 among all active vertices
    vertices[0] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);

    // decode e0

    // which of the edges incident to v0 is e0
    candidate_count = 0;
    // how many potential candidates are around that vertex?
    for (i = 0; i < vertices[0]->list_size; i++)
    {
      if (vertices[0]->list[i]->target == vertices[0])
      {
        edges[0] = vertices[0]->list[i];
        candidate_count++;
      }
    }
    // if there is more than one we need to decode which
    if (candidate_count > 1)
    {
      candidate = rd_conn->decode(candidate_count);
      PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
      for (i = 0; i < vertices[0]->list_size; i++)
      {
        if (vertices[0]->list[i]->target == vertices[0])
        {
          if (candidate)
          {
            candidate--;
          }
          else
          {
            edges[0] = vertices[0]->list[i];
            break;
          }
        }
      }
    }
    vertices[1] = edges[0]->origin;

    if (op == SMC_ADD)
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_add++;
#endif
      // allocate new vertex
      vertices[2] = allocVertex();
      // give it its index
      vertices[2]->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(vertices[0]->v, edges[0]->across, vertices[1]->v, vertices[2]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[2]);
      new_vertices[have_new] = 2;
      have_new++;

      // two edges are not known
      edges[1] = 0;
      edges[2] = 0;
    }
    else if (op == SMC_JOIN)
    {
#ifdef PRINT_CONTROL_OUTPUT
      op_join++;
#endif
       // decode v2's relative index 
      dv_index = rd_conn_index->decode(dv->size());
      // get v1 among all active vertices
      vertices[2] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      // two edges are not known
      edges[1] = 0;
      edges[2] = 0;
    }
    else // fill or end
    {
      // decode edge e2
      
      // which of the edges incident to v0 is e2
      candidate_count = 0;
      // how many potential candidates are around that vertex?
      for (i = 0; i < vertices[0]->list_size; i++)
      {
        if (vertices[0]->list[i]->origin == vertices[0])
        {
          edges[2] = vertices[0]->list[i];
          candidate_count++;
        }
      }
      // if there is more than one we need to decode which
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
        PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
        for (i = 0; i < vertices[0]->list_size; i++)
        {
          if (vertices[0]->list[i]->origin == vertices[0])
          {
            if (candidate)
            {
              candidate--;
            }
            else
            {
              edges[2] = vertices[0]->list[i];
              break;
            }
          }
        }
      }
      vertices[2] = edges[2]->target;

      // search for a potential e1 edge (e.g. it was an end)
      edges[1] = 0;
      for (i = 0; i < vertices[2]->list_size; i++)
      {
        if (vertices[2]->list[i]->target == vertices[1])
        {
          edges[1] = vertices[2]->list[i];
          break;
        }
      }
#ifdef PRINT_CONTROL_OUTPUT
      if (edges[1])
      {
        op_end++;
      }
      else
      {
        op_fill++;
      }
#endif
    }
  }

  // copy vertex and triangle data into API

  for (i = 0; i < 3; i++)
  {
    t_idx[i] = vertices[i]->index;
    t_pos_f[i] = vertices[i]->v;
  }

  // increment vertex use_counts, create edges, and update edge degrees

  for (i = 0; i< 3; i++)
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

  // write finalization info, finalize data structures

  SMvertex* vertex;
  SMedge* edge;
  int degree_one;
  int use_count;

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

    if (rd_conn_final->decode(rmFinalized[degree_one][use_count]))
    {
      t_final[i] = true;
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
      finalized_vertices[have_finalized] = i;
      have_finalized++;
      deallocVertex(vertex);
    }
    else
    {
      t_final[i] = false;
    }
  }
  return true;
}

bool SMreader_smd::decompress_triangle_traversal()
{
  int i,j;
  int op;
  int candidate;
  int candidate_count;

  SMedge* edges[3];
  SMvertex* vertices[3];

  edges[0] = (SMedge*)traversal_queue->getAndRemoveFirstElement();

  vertices[0] = edges[0]->target;
  vertices[1] = edges[0]->origin;

  edges[1] = 0;
  edges[2] = 0;

  i = (vertices[0]->use_count < MAX_USE_COUNT_OP)?vertices[0]->use_count:MAX_USE_COUNT_OP-1;
  j = (vertices[1]->use_count < MAX_USE_COUNT_OP)?vertices[1]->use_count:MAX_USE_COUNT_OP-1;
  op = rd_conn_op->decode(rmTraversalOp[i][j]); 

  PRINT_DEBUG_OUTPUT(stderr, "t%d:%d:%d ", i,j,op);

  if (op == SMC_ADD)
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_add++;
#endif
    // allocate new vertex
    vertices[2] = allocVertex();
    // give it its index
    vertices[2]->index = v_count + have_new;
    // decode v2's position
    decompressVertexPosition(vertices[0]->v, edges[0]->across, vertices[1]->v, vertices[2]->v);
    // insert it into the indexable data structure
    dv->addElement(vertices[2]);
    new_vertices[have_new] = 2;
    have_new++;
  }
  else if (op == SMC_FILL_END)
  {
    if (vertices[1]->use_count >= vertices[0]->use_count)
    {
      if (rd_conn_rl->decode(rmRight) == 0) // confirm
      {
        PRINT_DEBUG_OUTPUT(stderr, "r0 ");
#ifdef PRINT_CONTROL_OUTPUT
        right_confirm++;
#endif
        i = 1;
      }
      else // correct
      {
        PRINT_DEBUG_OUTPUT(stderr, "r1 ");
#ifdef PRINT_CONTROL_OUTPUT
        right_correct++;
#endif
        i = 0;
      }
    }
    else
    {
      if (rd_conn_rl->decode(rmLeft) == 0) // confirm
      {
        PRINT_DEBUG_OUTPUT(stderr, "l0 ");
#ifdef PRINT_CONTROL_OUTPUT
        left_confirm++;
#endif
        i = 0;
      }
      else // correct
      {
        PRINT_DEBUG_OUTPUT(stderr, "l1 ");
#ifdef PRINT_CONTROL_OUTPUT
        left_correct++;
#endif
        i = 1;
      }
    }
    if (i == 1)
    {
      // decode which is the edge e1 around v1 that gets us to v2
      candidate_count = 0;

      // how many edges are around v1 that could potentially be e1?
      for (i = 0; i < vertices[1]->list_size; i++)
      {
        if (vertices[1]->list[i]->target == vertices[1])
        {
          edges[1] = vertices[1]->list[i];
          candidate_count++;
        }
      }
      // only if there is more than one we need to decode it
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
        PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
        for (i = 0; i < vertices[1]->list_size; i++)
        {
          if (vertices[1]->list[i]->target == vertices[1])
          {
            if (candidate)
            {
              candidate--;
            }
            else
            {
              edges[1] = vertices[1]->list[i];
              break;
            }
          }
        }
      }
      vertices[2] = edges[1]->origin;
      // search for a potential e2 edge (e.g. it was an end)
      for (i = 0; i < vertices[2]->list_size; i++)
      {
        if (vertices[2]->list[i]->origin == vertices[0])
        {
          edges[2] = vertices[2]->list[i];
          break;
        }
      }
#ifdef PRINT_CONTROL_OUTPUT
      if (edges[2])
      {
        op_end++;
      }
      else
      {
        op_fill++;
      }
#endif
    }
    else if (i == 0)
    {
      candidate_count = 0;
      
      // how many edges are around v0 that could potentially be e2?
      for (i = 0; i < vertices[0]->list_size; i++)
      {
        if (vertices[0]->list[i]->origin == vertices[0])
        {
          edges[2] = vertices[0]->list[i];
          candidate_count++;
        }
      }
      // only if there is more than one we need to decode it
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
        PRINT_DEBUG_OUTPUT(stderr, "c%d:%d ",candidate_count,candidate);
        for (i = 0; i < vertices[0]->list_size; i++)
        {
          if (vertices[0]->list[i]->origin == vertices[0])
          {
            if (candidate)
            {
              candidate--;
            }
            else
            {
              edges[2] = vertices[0]->list[i];
              break;
            }
          }
        }
      }
      vertices[2] = edges[2]->target;
      // search for a potential e1 edge (e.g. it was an end)
      for (i = 0; i < vertices[2]->list_size; i++)
      {
        if (vertices[2]->list[i]->target == vertices[1])
        {
          edges[1] = vertices[2]->list[i];
          break;
        }
      }
#ifdef PRINT_CONTROL_OUTPUT
      if (edges[1])
      {
        op_end++;
      }
      else
      {
        op_fill++;
      }
#endif
    }
  }
  else if (op == SMC_JOIN)
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_join++;
#endif
    // decode v2's relative index
    i = rd_conn_index->decode(dv->size());        
    // get v2 among all active vertices
    vertices[2] = (SMvertex*)dv->getElementWithRelativeIndex(i);
  }
  else if (op == SMC_BORDER)
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_border++;
#endif
    removeEdgeFromVertex(edges[0],vertices[0]);
    removeEdgeFromVertex(edges[0],vertices[1]);
    deallocEdge(edges[0]);
    if (traversal_queue->elements())
    {
      return decompress_triangle_traversal();
    }
    return false;
  }
  else if (op == SMC_SKIP)
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_skip++;
#endif
    ((int*)(edges[0]))[0] = -1; // mark edge as no longer being in the queue
    return false;
  }

  // copy vertex and triangle data into API

  for (i = 0; i < 3; i++)
  {
    t_idx[i] = vertices[i]->index;
    t_pos_f[i] = vertices[i]->v;
  }

  // edges[0]

  removeEdgeFromVertex(edges[0],vertices[0]);
  removeEdgeFromVertex(edges[0],vertices[1]);
  deallocEdge(edges[0]);

  // edges[1]

  if (edges[1])
  {
    if (((int*)edges[1])[0] >= 0)
    {
      traversal_queue->removeElement(edges[1]);
    }
    removeEdgeFromVertex(edges[1],vertices[1]);
    removeEdgeFromVertex(edges[1],vertices[2]);
    deallocEdge(edges[1]);
  }
  else
  { 
    edges[1] = allocEdge(vertices[0]->v);  // the vertex coordinates are used for parallelogram predictions
    edges[1]->origin = vertices[1];
    edges[1]->target = vertices[2];
    addEdgeToVertex(edges[1],vertices[1]);
    addEdgeToVertex(edges[1],vertices[2]);
    traversal_queue->addElement(edges[1]);
  }

  // edges[2]

  if (edges[2])
  {
    if (((int*)edges[2])[0] >= 0)
    {
      traversal_queue->removeElement(edges[2]);
    }
    removeEdgeFromVertex(edges[2],vertices[2]);
    removeEdgeFromVertex(edges[2],vertices[0]);
    deallocEdge(edges[2]);
  }
  else
  { 
    edges[2] = allocEdge(vertices[1]->v);  // the vertex coordinates are used for parallelogram predictions
    edges[2]->origin = vertices[2];
    edges[2]->target = vertices[0];
    addEdgeToVertex(edges[2],vertices[2]);
    addEdgeToVertex(edges[2],vertices[0]);
    traversal_queue->addElement(edges[2]);
  }

  // increment vertex use_counts, write finalization info and finalize data structures

  SMvertex* vertex;
  SMedge* edge;
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
    
    if (rd_conn_final->decode(rmFinalized[degree_one][use_count])) // does this triangle finalize this vertex?
    {
      t_final[i] = true;
      dv->removeElement(vertex);
      for (j = 0; j < vertex->list_size; j++)
      {
        edge = vertex->list[j];
        if (vertex == edge->origin)
        {
          removeEdgeFromVertex(edge, edge->target);
          if (edge->target->list_size == 0)
          {
            // if not finalized it should be put into traversal queue
            // this check also needs to occur elsewhere
          }
        }
        else
        {
          removeEdgeFromVertex(edge, edge->origin);
          if (edge->origin->list_size == 0)
          {
            // if not finalized it should be put into traversal queue
            // this check also needs to occur elsewhere
          }
        }
        if (((int*)edge)[0] >= 0) traversal_queue->removeElement(edge);
        deallocEdge(edge);
      }
      finalized_vertices[have_finalized] = i;
      have_finalized++;
      deallocVertex(vertex);
    }
    else
    {
      t_final[i] = false;
    }
  }
  return true;
}

bool SMreader_smd::decompress_triangle()
{
  if (traversal_queue->elements())
  {
    if (next_waiting > 0)
    {
      if (decompress_triangle_traversal())
      {
        next_waiting--;
        have_triangle = 1;
        return true;
      }
    }
  }
    
  if (decompress_triangle_waiting())
  {
    next_waiting = (100 < traversal_queue->elements() ? traversal_queue->elements() : 100);
    have_triangle = 1;
    return true;
  }
  return false;
}

SMevent SMreader_smd::read_element()
{
  if (have_triangle == 0)
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
      have_triangle = -1;
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

SMevent SMreader_smd::read_event()
{
  if (have_triangle == 0 && have_finalized == 0)
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
      have_triangle = -1;
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
  else if (have_triangle > 0)
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
  else if (have_triangle == -1)
  {
    return SM_EOF;
  }
  return SM_ERROR;
}

SMreader_smd::SMreader_smd()
{
  // init of SMreader_smd
  nbits = -1;
  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;
}

SMreader_smd::~SMreader_smd()
{
  // clean-up for SMreader_smd
}
