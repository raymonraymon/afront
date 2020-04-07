/*
===============================================================================

  FILE:  smwriteopener.h
  
  CONTENTS:
  
    Takes a file_name and/or a format description as input and opens streaming
    (and non-streaming) meshes in various formats for write access. This avoids
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

    23 May 2007 -- created after picking house meeting food at house of curries

===============================================================================
*/
#ifndef SMWRITEOPENER_H
#define SMWRITEOPENER_H

#include "smwriter.h"

#include <stdio.h>

class SMwriteOpener
{
public:

  void set_file_format(const char* file_format);
  void set_file_name(const char* file_name);
  void set_compression_bits(int compression_bits);
  void set_maximum_delay(int maximum_delay);
  void set_dry_run(bool dry_run);

  SMwriter* open();

  SMwriteOpener();
  ~SMwriteOpener();

  char* file_format;
  char* file_name;
private:
  int compression_bits;
  int maximum_delay;
  bool dry_run;
};

#endif
