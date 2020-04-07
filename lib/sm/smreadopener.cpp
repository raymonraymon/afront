/*
===============================================================================

  FILE:  smreadopener.cpp
  
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
#include "filewrapper.h"
#include "smreadopener.h"

// streaming mesh readers
#include "smreader_sma.h"
#include "smreader_smb.h"
#include "smreader_smc.h"
#include "smreader_smd.h"
// only supported for legacy reasons
#include "smreader_smc_old.h"

// non-streaming mesh readers
#include "smreader_ply.h"
#include "smreader_off.h"
#include "smreader_vtk.h"
#include "smreader_jrs.h"

#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
extern "C" FILE* fopenGzipped(const char* filename, const char* mode);
#endif

void SMreadOpener::set_file_name(const char* file_name)
{
  if (this->file_name) free(this->file_name);
  this->file_name = strdup(file_name);
}

void SMreadOpener::set_file_format(const char* file_format)
{
  if (this->file_format) free(this->file_format);
  this->file_format = strdup(file_format);
}

SMreader* SMreadOpener::open()
{
  // if the file format is not given explicitely we try to guess it
  if (file_format == 0)
  {
    if (file_name == 0)
    {
      fprintf(stderr, "ERROR: neither a file name nor a file format were specified\n");
      return 0;
    }
    else if (strstr(file_name, ".sma") || strstr(file_name, ".obj"))
      file_format = strdup("sma");
    else if (strstr(file_name, ".smb"))
      file_format = strdup("smb");
    else if (strstr(file_name, ".smc_old"))
      file_format = strdup("smc_old");
    else if (strstr(file_name, ".smc"))
      file_format = strdup("smc");
    else if (strstr(file_name, ".smd"))
      file_format = strdup("smd");
    else if (strstr(file_name, ".ply"))
      file_format = strdup("ply");
    else if (strstr(file_name, ".off"))
      file_format = strdup("off");
    else if (strstr(file_name, ".vtk"))
      file_format = strdup("vtk");
    else if (strstr(file_name, ".node"))
      file_format = strdup("jrs");
    else
    {
      fprintf(stderr, "ERROR: could not guess file format from file name '%s'\n", file_name);
      return 0;
    }
  }

  // depending on the file format we open the file in binary or text mode
  bool is_binary = true;
  if (strcmp(file_format, "sma") == 0 || strcmp(file_format, "jrs") == 0 || strcmp(file_format, "off") == 0 || strcmp(file_format, "vtk") == 0)
    is_binary = false;

  FILE* file;

  // if we have no file name we read from stdin, otherwise we open the file
  if (file_name == 0)
  {
	  file = myfopen(stdin, is_binary ? "wb" : "w");
#ifdef _WIN32
    if (is_binary)
      if(_setmode( _fileno( stdin ), _O_BINARY ) == -1 )
      {
        fprintf(stderr, "ERROR: cannot set stdin to binary (untranslated) mode\n");
      }
    else
      if(_setmode( _fileno( stdin ), _O_TEXT ) == -1 )
      {
        fprintf(stderr, "ERROR: cannot set stdin to text (translated) mode\n");
      }
#endif
  }
  else
  {
    if (strstr(file_name, ".gz"))
    {
#ifdef _WIN32
      file = fopenGzipped(file_name, (is_binary ? "rb" : "r"));
#else
	  file = myfopen(file_name, (is_binary ? "rb" : "r"));
	  //      fprintf(stderr,"ERROR: cannot open gzipped file '%s'\n",file_name);
	  //      return 0;
#endif
    }
    else
    {
      file = myfopen(file_name, (is_binary ? "rb" : "r"));
    }
    if (file == 0)
    {
      fprintf(stderr,"ERROR: cannot open file '%s'\n",file_name);
      return 0;
    }
  }

  // now we open the respective reader

  if (strcmp(file_format, "sma") == 0)
  {
    SMreader_sma* smreader_sma = new SMreader_sma();
    if (smreader_sma->open(file)) return smreader_sma;
    delete smreader_sma;
  }
  else if (strcmp(file_format, "smb") == 0)
  {
    SMreader_smb* smreader_smb = new SMreader_smb();
    if (smreader_smb->open(file)) return smreader_smb;
    delete smreader_smb;
  }
  else if (strcmp(file_format, "smc") == 0)
  {
    SMreader_smc* smreader_smc = new SMreader_smc();
    if (smreader_smc->open(file)) return smreader_smc;
    delete smreader_smc;
  }
  else if (strcmp(file_format, "smd") == 0)
  {
    SMreader_smd* smreader_smd = new SMreader_smd();
    if (smreader_smd->open(file)) return smreader_smd;
    delete smreader_smd;
  }
  else if (strcmp(file_format, "smc_old") == 0)
  {
    SMreader_smc_old* smreader_smc_old = new SMreader_smc_old();
    if (smreader_smc_old->open(file)) return smreader_smc_old;
    delete smreader_smc_old;
  }
  else if (strcmp(file_format, "ply") == 0)
  {
    SMreader_ply* smreader_ply = new SMreader_ply();
    if (smreader_ply->open(file)) return smreader_ply;
    delete smreader_ply;
  }
  else if (strcmp(file_format, "off") == 0)
  {
    SMreader_off* smreader_off = new SMreader_off();
    if (smreader_off->open(file)) return smreader_off;
    delete smreader_off;
  }
  else if (strcmp(file_format, "vtk") == 0)
  {
    SMreader_vtk* smreader_vtk = new SMreader_vtk();
    if (smreader_vtk->open(file)) return smreader_vtk;
    delete smreader_vtk;
  }
  else if (strcmp(file_format, "jrs") == 0)
  {
    FILE* file_other;
    char* file_name_other = strdup(file_name);
    int len = strlen(file_name);
    do len--; while (len > 0 && file_name_other[len] != '.');
    if (strstr(file_name, ".gz"))
    {
      do len--; while (len > 0 && file_name_other[len] != '.');
      file_name_other[len+1] = 'e';
      file_name_other[len+2] = 'l';
      file_name_other[len+3] = 'e';
      file_name_other[len+4] = '.';
      file_name_other[len+5] = 'g';
      file_name_other[len+6] = 'z';
      file_name_other[len+7] = '\0';
#ifdef _WIN32
      file_other = fopenGzipped(file_name_other, "r");
#else      
      file_other = 0;
#endif
    }
    else
    {
      file_name_other[len+1] = 'e';
      file_name_other[len+2] = 'l';
      file_name_other[len+3] = 'e';
      file_name_other[len+4] = '\0';
      file_other = myfopen(file_name_other, "r");
    }
    if (file_other == 0)
    {
      fprintf(stderr,"ERROR: cannot open other file '%s'\n",file_name);
      return 0;
    }
    SMreader_jrs* smreader_jrs = new SMreader_jrs();
    if (smreader_jrs->open(file,file_other)) return smreader_jrs;
    delete smreader_jrs;
  }
  else
  {
    fprintf(stderr,"ERROR: cannot determine the reader for format '%s'\n", file_format);
   }
  return 0;
}

bool SMreadOpener::reopen(SMreader* smreader)
{
  if (file_name == 0)
  {
    fprintf(stderr, "ERROR: cannot reopen() without a file name\n");
    return false;
  }
  else if (file_format == 0)
  {
    fprintf(stderr, "ERROR: cannot reopen() without a successful open()\n");
    return false;
  }
  else if (smreader == 0)
  {
    fprintf(stderr, "ERROR: cannot reopen() when smreader is a null pointer\n");
    return false;
  }

  // close the reader
  smreader->close();

  // depending on the file format we open the file in binary or text mode
  bool is_binary = true;
  if (strcmp(file_format, "sma") == 0 || strcmp(file_format, "off") == 0 || strcmp(file_format, "vtk") == 0 || strcmp(file_format, "jrs") == 0)
    is_binary = false;

  FILE* file;

  // reopen the file
  if (strstr(file_name, ".gz"))
  {
#ifdef _WIN32
    file = fopenGzipped(file_name, (is_binary ? "rb" : "r"));
#else
    file = myfopen(file_name, (is_binary ? "rb" : "r"));
	//    fprintf(stderr,"ERROR: cannot reopen gzipped file '%s'\n",file_name);
	//    return false;
#endif
  }
  else
  {
    file = myfopen(file_name, (is_binary ? "rb" : "r"));
  }
  if (file == 0)
  {
    fprintf(stderr,"ERROR: cannot reopen file '%s'\n",file_name);
    return false;
  }

  // reopen the respective reader
  if (strcmp(file_format, "sma") == 0)
  {
    SMreader_sma* smreader_sma = (SMreader_sma*)smreader;
    return smreader_sma->open(file);
  }
  else if (strcmp(file_format, "smb") == 0)
  {
    SMreader_smb* smreader_smb = (SMreader_smb*)smreader;
    return smreader_smb->open(file);
  }
  else if (strcmp(file_format, "smc") == 0)
  {
    SMreader_smc* smreader_smc = (SMreader_smc*)smreader;
    return smreader_smc->open(file);
  }
  else if (strcmp(file_format, "smd") == 0)
  {
    SMreader_smd* smreader_smd = (SMreader_smd*)smreader;
    return smreader_smd->open(file);
  }
  else if (strcmp(file_format, "smc_old") == 0)
  {
    SMreader_smc_old* smreader_smc_old = (SMreader_smc_old*)smreader;
    return smreader_smc_old->open(file);
  }
  else if (strcmp(file_format, "ply") == 0)
  {
    SMreader_ply* smreader_ply = (SMreader_ply*)smreader;
    return smreader_ply->open(file);
  }
  else if (strcmp(file_format, "off") == 0)
  {
    SMreader_off* smreader_off = (SMreader_off*)smreader;
    return smreader_off->open(file);
  }
  else if (strcmp(file_format, "vtk") == 0)
  {
    SMreader_vtk* smreader_vtk = (SMreader_vtk*)smreader;
    return smreader_vtk->open(file);
  }
  else if (strcmp(file_format, "jrs") == 0)
  {
    FILE* file_other;
    char* file_name_other = strdup(file_name);
    int len = strlen(file_name);
    do len--; while (len > 0 && file_name_other[len] != '.');
    if (strstr(file_name, ".gz"))
    {
      do len--; while (len > 0 && file_name_other[len] != '.');
      file_name_other[len+1] = 'e';
      file_name_other[len+2] = 'l';
      file_name_other[len+3] = 'e';
      file_name_other[len+4] = '.';
      file_name_other[len+5] = 'g';
      file_name_other[len+6] = 'z';
      file_name_other[len+7] = '\0';
#ifdef _WIN32
      file_other = fopenGzipped(file_name_other, "r");
#else      
      fprintf(stderr,"ERROR: cannot reopen other gzipped file '%s'\n",file_name_other);
      return false;
#endif
    }
    else
    {
      file_name_other[len+1] = 'e';
      file_name_other[len+2] = 'l';
      file_name_other[len+3] = 'e';
      file_name_other[len+4] = '\0';
      file_other = myfopen(file_name_other, "r");
    }
    if (file_other == 0)
    {
      fprintf(stderr,"ERROR: cannot reopen other file '%s'\n",file_name);
      return false;
    }
    SMreader_jrs* smreader_jrs = (SMreader_jrs*)smreader;
    return smreader_jrs->open(file,file_other);
  }
  else
  {
    fprintf(stderr,"ERROR: cannot determine the reader for format '%s'\n", file_format);
  }
  return false;
}

SMreadOpener::SMreadOpener()
{
  file_name = 0;
  file_format = 0;
}

SMreadOpener::~SMreadOpener()
{
  if (file_format) free(file_format);
  if (file_name) free(file_name);
}
