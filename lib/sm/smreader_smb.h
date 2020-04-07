/*
===============================================================================

  FILE:  smreader_smb.h
  
  CONTENTS:
   
    Reads a Streaming Mesh from an efficient binary format. 
    
    Optionally the connectivity and the geometry may be compressed using a
    simple lossless format that preserves both, the exact ordering of the
    mesh elements and the rotation of the triangles. But this is not yet
    implemented.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2004  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    1 August 2004 -- initial version created outside at Weaver Street Market
  
===============================================================================
*/
#ifndef SMREADER_SMB_H
#define SMREADER_SMB_H

#include "smreader.h"

#include <stdio.h>
#include "filewrapper.h"

class SMreader_smb : public SMreader
{
public:

  // smreader interface function implementations

  SMevent read_element();
  SMevent read_event();

  void close(bool close_file=true); // virtual base function must be reimplemented

  // smreader_sma functions

  bool open(FILE* fp);

  SMreader_smb();
  ~SMreader_smb();

private:
  FILE* file;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  void read_header();
  void read_buffer();

  bool endian_swap;

  int element_number;
  int element_counter;
  unsigned int element_descriptor;
  int* element_buffer;
};

#endif
