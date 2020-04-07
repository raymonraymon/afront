/*
===============================================================================

  FILE:  smreadopener.h
  
  CONTENTS:
  
    Takes a file_name and/or a format description as input and opens streaming
    (and non-streaming) meshes in various formats for read access. This avoids to
    (a) having to reimplement opening different formats in new applications
    (b) having to add new formats to every existing application

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2007  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:

    14 June 2007 -- stripped of smreader->close() which destroys filter paradigm
    23 May 2007 -- created after picking up an IKEA shelf to turn into firewood

===============================================================================
*/
#ifndef SMREADOPENER_H
#define SMREADOPENER_H

#include "smreader.h"

#include <stdio.h>

class SMreadOpener
{
public:

  void set_file_name(const char* file_name);
  void set_file_format(const char* file_format);

  SMreader* open();
  bool reopen(SMreader* smreader);

  SMreadOpener();
  ~SMreadOpener();

  char* file_name;
  char* file_format;
};

#endif
