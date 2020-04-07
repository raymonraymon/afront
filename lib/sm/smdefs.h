/*
===============================================================================

  FILE:  smdefs.h
  
  CONTENTS:
  
    Common definitions for our Streaming Mesh classes
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2007  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    14 June 2007 -- created after rejecting the offer from NU Singapore
  
===============================================================================
*/
#ifndef SMDEFS_H
#define SMDEFS_H

// events 
typedef enum {
  SM_ERROR = -1,
  SM_EOF = 0,
  SM_VERTEX = 1,
  SM_TRIANGLE = 2,
  SM_FINALIZED = 3
} SMevent;

// order 
typedef enum {
  SM_NO_ORDER= 0,
  SM_PRE_ORDER = 1,
  SM_POST_ORDER = 2
} SMorder;

#endif
