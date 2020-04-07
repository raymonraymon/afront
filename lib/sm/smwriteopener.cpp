/*
===============================================================================

  FILE:  smwriteopener.cpp
  
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
#include "smwriteopener.h"

// streaming mesh writers
#include "smwriter_sma.h"
#include "smwriter_smb.h"
#include "smwriter_smc.h"
#include "smwriter_smd.h"

// non-streaming mesh writers
#include "smwriter_ply.h"
#include "smwriter_off.h"

// for buffered writing
#include "smwritebuffered.h"

#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#include <io.h>
#endif

#include <fcntl.h>

#include "filewrapper.h"



void SMwriteOpener::set_file_format(const char* file_format)
{
  if (this->file_format) free(this->file_format);
  this->file_format = strdup(file_format);
}

void SMwriteOpener::set_file_name(const char* file_name)
{
  if (this->file_name) free(this->file_name);
  this->file_name = strdup(file_name);
}

void SMwriteOpener::set_compression_bits(int compression_bits)
{
  this->compression_bits = compression_bits;
}

void SMwriteOpener::set_maximum_delay(int maximum_delay)
{
  this->maximum_delay = maximum_delay;
}

void SMwriteOpener::set_dry_run(bool dry_run)
{
  this->dry_run = dry_run;
}

SMwriter* SMwriteOpener::open()
{
  // if the file format is not given explicitely we try to guess it
  if (file_format == 0)
  {
    if (file_name == 0)
    {
      fprintf(stderr, "ERROR: neither a file name nor a file format were specified\n");
      return 0;
    }
    else if (strstr(file_name, ".sma"))
      file_format = strdup("sma");
    else if (strstr(file_name, ".smb"))
      file_format = strdup("smb");
    else if (strstr(file_name, ".smc"))
      file_format = strdup("smc");
    else if (strstr(file_name, ".smd"))
      file_format = strdup("smd");
    else if (strstr(file_name, ".ply"))
      file_format = strdup("ply");
    else if (strstr(file_name, ".off"))
      file_format = strdup("off");
    else
    {
      fprintf(stderr, "ERROR: could not guess file format from file name '%s'\n", file_name);
      return 0;
    }
  }

  // depending on the file format we open the file in binary or text mode
  bool is_binary = true;
  if (strcmp(file_format, "sma") == 0 || strcmp(file_format, "off") == 0)
    is_binary = false;

  FILE* file;

  // if we have no file name we write to stdout, otherwise we open the file
  if (file_name == 0)
  {
	  file = myfopen(stdout, is_binary ? "wb" : "w");
#ifdef _WIN32
    if (is_binary)
      if(_setmode( _fileno( stdout ), _O_BINARY ) == -1 )
      {
        fprintf(stderr, "ERROR: cannot set stdout to binary (untranslated) mode\n");
      }
    else
      if(_setmode( _fileno( stdout ), _O_TEXT ) == -1 )
      {
        fprintf(stderr, "ERROR: cannot set stdout to text (translated) mode\n");
      }
#endif
  }
  else
  {
    if (strstr(file_name, ".gz"))
    {
#ifdef _WIN32
      //file = fopenGzipped(file_name, (is_binary ? "wb" : "w"));
#else
	  file = myfopen(file_name, (is_binary ? "wb" : "w"));
	  //      fprintf(stderr,"ERROR: cannot open gzipped file '%s'\n",file_name);
	  //      return 0;
#endif
    }
    else
    {
      file = myfopen(file_name, (is_binary ? "wb" : "w"));
    }
    if (file == 0)
    {
      fprintf(stderr,"ERROR: cannot open file '%s'\n",file_name);
      return 0;
    }
  }

  // now we open the respective writer

  if (strcmp(file_format, "sma") == 0)
  {
    SMwriter_sma* smwriter_sma = new SMwriter_sma();
    if (smwriter_sma->open(file)) return smwriter_sma;
    delete smwriter_sma;
  }
  else if (strcmp(file_format, "smb") == 0)
  {
    SMwriter_smb* smwriter_smb = new SMwriter_smb();
    if (smwriter_smb->open(file)) return smwriter_smb;
    delete smwriter_smb;
  }
  else if (strcmp(file_format, "smc") == 0)
  {
    SMwriter_smc* smwriter_smc = new SMwriter_smc();
    if (smwriter_smc->open(dry_run ? 0 : file, compression_bits))
    {
      if (maximum_delay)
      {
        SMwriteBuffered* smwrite_buffered = new SMwriteBuffered();
        if (smwrite_buffered->open(smwriter_smc, maximum_delay)) return smwrite_buffered;
        delete smwrite_buffered;
      }
      else
      {
        return smwriter_smc;
      }
    }
    else
    {
      delete smwriter_smc;
    }
  }
  else if (strcmp(file_format, "smd") == 0)
  {
    SMwriter_smd* smwriter_smd = new SMwriter_smd();
    if (maximum_delay)
    {
      if (smwriter_smd->open(dry_run ? 0 : file, compression_bits, maximum_delay)) return smwriter_smd;
    }
    else
    {
      if (smwriter_smd->open(dry_run ? 0 : file, compression_bits)) return smwriter_smd;
    }
    delete smwriter_smd;
  }
  else if (strcmp(file_format, "ply") == 0)
  {
    SMwriter_ply* smwriter_ply = new SMwriter_ply();
    if (smwriter_ply->open(file)) return smwriter_ply;
    delete smwriter_ply;
  }
  else if (strcmp(file_format, "off") == 0)
  {
    SMwriter_off* smwriter_off = new SMwriter_off();
    if (smwriter_off->open(file)) return smwriter_off;
    delete smwriter_off;
  }
  else
  {
    fprintf(stderr,"ERROR: cannot determine the writer for format '%s'\n", file_format);
  }
  return 0;
}

SMwriteOpener::SMwriteOpener()
{
  file_format = 0;
  file_name = 0;
  compression_bits = 16;
  maximum_delay = 10000;
  dry_run = false;
}

SMwriteOpener::~SMwriteOpener()
{
  if (file_format) free(file_format);
  if (file_name) free(file_name);
}
