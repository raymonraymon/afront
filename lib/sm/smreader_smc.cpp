/*
===============================================================================

  FILE:  smreader_smc.cpp
  
  CONTENTS:
  
    see header file

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see header file
  
===============================================================================
*/
#include "smreader_smc.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "rangemodel.h"
#include "rangedecoder.h"

#include "dynamicvector.h"
#include "littlecache.h"

#include "floatcompressor.h"

#include "positionquantizer_new.h"
#include "integercompressor_new.h"

#include "vec3fv.h"
#include "vec3iv.h"

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
  SMvertex* buffer_next; // used for efficient memory management (and by the dynamic vector)
  float v[3];
  int index;
  int use_count;
  int list_size;
  int list_alloc;
  SMedge** list;
} SMvertex;

typedef struct SMedge
{
  SMedge* buffer_next; // used for efficient memory management
  SMvertex* origin;
  SMvertex* target;
  float across[3];
} SMedge;

static DynamicVector* dv;

static LittleCache* lc;

static FloatCompressor* fc[3];

static PositionQuantizerNew* pq;
static IntegerCompressorNew* ic[3];

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

// we need to handle both versions SME and SME_NON_FINALIZED_EOF
static int version;

// what was the last operation
static int last_op;

// codes next operation
static RangeModel** rmOp;

// codes non-manifoldness of start operations
static RangeModel* rmS_Old;

// codes cache hits for start operations
static RangeModel** rmS_Cache;

// codes add/join operations
static RangeModel** rmAJ_Cache;

// codes fill/end operations
static RangeModel** rmFE_Cache;

// codes vertex finalization
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
  rmOp = (RangeModel**)malloc(sizeof(RangeModel*)*5);
  rmOp[0] = new RangeModel(4,0,compress); // after SMC_START
  rmOp[1] = new RangeModel(4,0,compress); // after SMC_ADD
  rmOp[2] = new RangeModel(4,0,compress); // after SMC_JOIN
  rmOp[3] = new RangeModel(4,0,compress); // after SMC_FILL
  rmOp[4] = new RangeModel(4,0,compress); // after SMC_END

  // range tables for start (S) operation
  if (version == SM_VERSION_SME_NON_FINALIZED_EOF)
  {
    rmS_Old = new RangeModel(5,0,compress); // number of old (e.g. non-manifold) vertices in a start operation (4 = NON_FINALIZED_EOF)
  }
  else
  {
    rmS_Old = new RangeModel(4,0,compress); // number of old (e.g. non-manifold) vertices in a start operation
  }

  // range tables for start (S) operations

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

// statistics

#ifdef PRINT_CONTROL_OUTPUT

static int op_start = 0;
static int op_add_join = 0;
static int op_fill_end = 0;

static int add_miss = 0;
static int add_hit[] = {0,0,0,0,0,0};

static int fill_miss = 0;
static int fill_hit[] = {0,0,0,0,0,0,0,0,0};

static int start_non_manifold = 0;
static int add_non_manifold = 0;

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

bool SMreader_smc::open(FILE* file)
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

  version = fgetc(file);
  // read version
  if (version != SM_VERSION_SME && version != SM_VERSION_SME_NON_FINALIZED_EOF)
  {
    fprintf(stderr,"ERROR: wrong SMreader (need %d but this is SMreader_smc (%d or %d)\n",version,SM_VERSION_SME,SM_VERSION_SME_NON_FINALIZED_EOF);
    if (version == EOF) fprintf(stderr,"       check ... maybe the file (or the stream) is empty\n");
    exit(0);
  }

  this->file = file;

  initVertexBuffer(1024);
  initEdgeBuffer(1024);

  dv = new DynamicVector();
  lc = new LittleCache();

  initDecoder(file);
  initModels(0);

  last_op = 0;

  // read precision
  nbits = rd_conn->decode(25);

  v_count = 0;
  f_count = 0;

  read_header();

  return true;
}

void SMreader_smc::close(bool close_file)
{
  // close of SMreader_smc
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

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"edge_buffer_size %d edge_buffer_maxsize %d\n",edge_buffer_size,edge_buffer_maxsize);
  fprintf(stderr,"vertex_buffer_size %d vertex_buffer_maxsize %d\n",vertex_buffer_size,vertex_buffer_maxsize);
  fprintf(stderr,"op_start %d op_add_join %d op_fill_end %d\n",op_start,op_add_join,op_fill_end);
  fprintf(stderr,"start_non_manifold %d add_non_manifold %d\n",start_non_manifold,add_non_manifold);
  fprintf(stderr,"add_miss %d add_hit %d (%d %d %d %d %d %d)\n",add_miss,add_hit[0]+add_hit[1]+add_hit[2]+add_hit[3]+add_hit[4]+add_hit[5],add_hit[0],add_hit[1],add_hit[2],add_hit[3],add_hit[4],add_hit[5]);
  fprintf(stderr,"fill_miss %d fill_hit %d (%d %d %d %d %d %d %d %d %d)\n",fill_miss,fill_hit[0]+fill_hit[1]+fill_hit[2]+fill_hit[3]+fill_hit[4]+fill_hit[5]+fill_hit[6]+fill_hit[7]+fill_hit[8],fill_hit[0],fill_hit[1],fill_hit[2],fill_hit[3],fill_hit[4],fill_hit[5],fill_hit[6],fill_hit[7],fill_hit[8]);
#endif 

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

void SMreader_smc::read_header()
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

int SMreader_smc::decompress_triangle()
{
  int i,j;
  int op;
  int lc_pos;
  int dv_index;
  int degree_one;
  int use_count;
  int candidate;
  int candidate_count;

  SMvertex* vertices[3];
  SMedge* edges[3];

  if (edge_buffer_size)
  {
    op = rd_conn_op->decode(rmOp[last_op]);
  }
  else
  {
    if (rd_conn_op->decode(rmDone))
    {
      return 0;
    }
    else
    {
      op = SMC_START;
    }
  }

  if (op == SMC_ADD || op == SMC_JOIN)
  {
    // decode the index of v0 (-> improve with cache of with prediction)
    lc_pos = rd_conn_cache->decode(rmAJ_Cache[last_op]);
    if (lc_pos == 6)
    {
      lc_pos = -1;
#ifdef PRINT_CONTROL_OUTPUT
      add_miss++;
#endif
      dv_index = rd_conn_index->decode(dv->size());        
      vertices[0] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
    }
    else if (lc_pos < 3)
    {
#ifdef PRINT_CONTROL_OUTPUT
      add_hit[lc_pos]++;
#endif
      vertices[0] = (SMvertex*)lc->get(lc_pos);
    }
    else
    {
#ifdef PRINT_CONTROL_OUTPUT
      add_hit[lc_pos]++;
#endif
      vertices[1] = (SMvertex*)lc->get(lc_pos-3);
    }

    if (lc_pos < 3)
    {
      // decode which of the edges incident to v0 is e0
      candidate_count = 0;

      // how many potential candidate edges are around that vertex?
      for (i = 0; i < vertices[0]->list_size; i++)
      {
        if (vertices[0]->list[i]->target == vertices[0])
        {
          edges[0] = vertices[0]->list[i];
          candidate_count++;
        }
      }
      // only if there is more than one we need to decode it
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
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
    }
    else
    {
      // decode which of the edges incident to v1 is e0
      candidate_count = 0;

      // how many potential candidate edges are around that vertex?
      for (i = 0; i < vertices[1]->list_size; i++)
      {
        if (vertices[1]->list[i]->origin == vertices[1])
        {
          edges[0] = vertices[1]->list[i];
          candidate_count++;
        }
      }
      // only if there is more than one we need to decode it
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
        for (i = 0; i < vertices[1]->list_size; i++)
        {
          if (vertices[1]->list[i]->origin == vertices[1])
          {
            if (candidate)
            {
              candidate--;
            }
            else
            {
              edges[0] = vertices[1]->list[i];
              break;
            }
          }
        }
      }
      vertices[0] = edges[0]->target;
    }

    // encode v2
    if (op == SMC_JOIN)
    {
      // an old vertex
      dv_index = rd_conn_index->decode(dv->size());        
      vertices[2] = (SMvertex*) dv->getElementWithRelativeIndex(dv_index);
#ifdef PRINT_CONTROL_OUTPUT
      add_non_manifold++;
#endif
    }
    else
    {
      // a new vertex
      vertices[2] = allocVertex();
      // give it its index
      vertices[2]->index = v_count + have_new;
      // decode its position
      decompressVertexPosition(vertices[0]->v, edges[0]->across, vertices[1]->v, vertices[2]->v);
      // insert it into the indexable data structure
      dv->addElement(vertices[2]);
      new_vertices[have_new] = 2;
      have_new++;
    }
    edges[1] = 0;
    edges[2] = 0;
  }
  else if (op == SMC_FILL_END)
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_fill_end++;
#endif

    // decode which active vertex we use for this triangle
    lc_pos = rd_conn_cache->decode(rmFE_Cache[last_op]);
    if (lc_pos == 9)
    {
      lc_pos = -1;
#ifdef PRINT_CONTROL_OUTPUT
      fill_miss++;
#endif
      dv_index = rd_conn_index->decode(dv->size());
      vertices[0] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
    }
    else if (lc_pos < 3)
    {
#ifdef PRINT_CONTROL_OUTPUT
      fill_hit[lc_pos]++;
#endif
      vertices[0] = (SMvertex*)lc->get(lc_pos);
    }
    else if (lc_pos < 6)
    {
#ifdef PRINT_CONTROL_OUTPUT
      fill_hit[lc_pos]++;
#endif
      vertices[1] = (SMvertex*)lc->get(lc_pos-3);
    }
    else
    {
#ifdef PRINT_CONTROL_OUTPUT
      fill_hit[lc_pos]++;
#endif
      vertices[2] = (SMvertex*)lc->get(lc_pos-6);
    }

    if (lc_pos < 3) // we have v0
    {
      // decode edge e0 (e.g. the edge from v0 to v1)
      candidate_count = 0;
      // how many edges are around v0 that could potentially be e0?
      for (i = 0; i < vertices[0]->list_size; i++)
      {
        if (vertices[0]->list[i]->target == vertices[0])
        {
          edges[0] = vertices[0]->list[i];
          candidate_count++;
        }
      }
      // only if there is more than one we need to decode it
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
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

      // decode edge e2 (e.g. the edge from v2 to v0)
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
    }
    else if (lc_pos < 6) // we have v1
    {
      // decode edge e0 (e.g. the edge from v0 to v1)
      candidate_count = 0;
      // how many edges are around v1 that could potentially be e0?
      for (i = 0; i < vertices[1]->list_size; i++)
      {
        if (vertices[1]->list[i]->origin == vertices[1])
        {
          edges[0] = vertices[1]->list[i];
          candidate_count++;
        }
      }
      // only if there is more than one we need to decode it
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
        for (i = 0; i < vertices[1]->list_size; i++)
        {
          if (vertices[1]->list[i]->origin == vertices[1])
          {
            if (candidate)
            {
              candidate--;
            }
            else
            {
              edges[0] = vertices[1]->list[i];
              break;
            }
          }
        }
      }
      vertices[0] = edges[0]->target;

      // decode edge e2 (e.g. the edge from v2 to v0)
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
    }
    else // we have v2
    {
      // decode edge e2 (e.g. the edge from v2 to v0)
      candidate_count = 0;
      // how many edges are around v2 that could potentially be e2?
      for (i = 0; i < vertices[2]->list_size; i++)
      {
        if (vertices[2]->list[i]->target == vertices[2])
        {
          edges[2] = vertices[2]->list[i];
          candidate_count++;
        }
      }
      // only if there is more than one we need to decode it
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
        for (i = 0; i < vertices[2]->list_size; i++)
        {
          if (vertices[2]->list[i]->target == vertices[2])
          {
            if (candidate)
            {
              candidate--;
            }
            else
            {
              edges[2] = vertices[2]->list[i];
              break;
            }
          }
        }
      }
      vertices[0] = edges[2]->origin;

      // decode edge e0 (e.g. the edge from v0 to v1)
      candidate_count = 0;
      // how many edges are around v0 that could potentially be e0?
      for (i = 0; i < vertices[0]->list_size; i++)
      {
        if (vertices[0]->list[i]->target == vertices[0])
        {
          edges[0] = vertices[0]->list[i];
          candidate_count++;
        }
      }
      // only if there is more than one we need to decode it
      if (candidate_count > 1)
      {
        candidate = rd_conn->decode(candidate_count);
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
    }

    // is it a FILL or an END operation

    edges[1] = 0;
    op = SMC_FILL;
    for (i = 0; i < vertices[2]->list_size; i++)
    {
      if (vertices[2]->list[i]->target == vertices[1])
      {
        edges[1] = vertices[2]->list[i];
        op = SMC_END;
        break;
      }
    }
  }
  else if (op == SMC_START)
  {
#ifdef PRINT_CONTROL_OUTPUT
    op_start++;
#endif

    // how many non-manifold vertices?
    use_count = rd_conn_op->decode(rmS_Old);

    if (use_count == 4)
    {
      return 0;
    }

    // decode vertex v0

    if (use_count)
    {
      // an old vertex
      lc_pos = rd_conn_cache->decode(rmS_Cache[0]);
      if (lc_pos == 3)
      {
        dv_index = rd_conn_index->decode(dv->size());        
        vertices[0] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      }
      else
      {
        vertices[0] = (SMvertex*)lc->get(lc_pos);
      }
#ifdef PRINT_CONTROL_OUTPUT
      start_non_manifold++;
#endif
      use_count--;
    }
    else
    {
      // a new vertex
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

    if (use_count)
    {
      // an old vertex
      lc_pos = rd_conn_cache->decode(rmS_Cache[1]);
      if (lc_pos == 3)
      {
        dv_index = rd_conn_index->decode(dv->size());        
        vertices[1] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      }
      else
      {
        vertices[1] = (SMvertex*)lc->get(lc_pos);
      }
#ifdef PRINT_CONTROL_OUTPUT
      start_non_manifold++;
#endif
      use_count--;
    }
    else
    {
      // a new vertex
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

    if (use_count)
    {
      // an old vertex
      lc_pos = rd_conn_cache->decode(rmS_Cache[2]);
      if (lc_pos == 3)
      {
        dv_index = rd_conn_index->decode(dv->size());        
        vertices[2] = (SMvertex*)dv->getElementWithRelativeIndex(dv_index);
      }
      else
      {
        vertices[2] = (SMvertex*)lc->get(lc_pos);
      }
#ifdef PRINT_CONTROL_OUTPUT
      start_non_manifold++;
#endif
      use_count--;
    }
    else
    {
      // a new vertex
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
    edges[0] = 0;
    edges[1] = 0;
    edges[2] = 0;
  }

  lc->put(vertices[0], vertices[1], vertices[2]);

  last_op = op;

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
      addEdgeToVertices(edges[i]);
    }
    else
    {
      removeEdgeFromVertices(edges[i]);
      deallocEdge(edges[i]);
    }
  }

  // read finalization info and finalize vertices

  for (i = 0; i < 3; i++)
  {
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
  }
  have_triangle = 1;
  return 1;
}

SMevent SMreader_smc::read_element()
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

SMevent SMreader_smc::read_event()
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

SMreader_smc::SMreader_smc()
{
  // init of SMreader_smc
  nbits = -1;
  have_new = 0; next_new = 0;
  have_triangle = 0;
  have_finalized = 0; next_finalized = 0;
}

SMreader_smc::~SMreader_smc()
{
  // clean-up for SMwriter_smc
}
