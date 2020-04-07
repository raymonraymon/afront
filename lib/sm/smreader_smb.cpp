/*
===============================================================================

  FILE:  smreader_smb.cpp
  
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
#include "smreader_smb.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "vec3fv.h"
#include "vec3iv.h"

#define SM_VERSION 0 // this is SMB

bool SMreader_smb::open(FILE* file)
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

  this->file = file;

  int input = fgetc(file);
  // read version
  if (input != SM_VERSION)
  {
    fprintf(stderr,"ERROR: wrong SMreader (need %d but this is SMreader_smb %d)\n",input,SM_VERSION);
    if (input == EOF) fprintf(stderr,"       check ... maybe the file (or the stream) is empty\n");
    exit(0);
  }

  read_header();
  read_buffer();

  if (element_descriptor & 1)
  {
    post_order = false;
  }
  else
  {
    post_order = true;
  }

  v_count = 0;
  f_count = 0;

  return true;
}

void SMreader_smb::close(bool close_file)
{
  // close of SMreader_smb
	if (close_file && file /*&& file != stdin*/) fclose(file);
  file = 0;
  have_finalized = 0; next_finalized = 0;

  element_number = 0;
  element_counter = 0;

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

SMevent SMreader_smb::read_element()
{
  if (element_counter < element_number)
  {
    have_finalized = next_finalized = 0;
    if (element_descriptor & 1) // next element is a vertex
    {
      v_idx = v_count;
      if (endian_swap) VecCopy3fv_swap_endian(v_pos_f, (float*)(&element_buffer[element_counter*3]));
      else VecCopy3fv(v_pos_f, (float*)(&element_buffer[element_counter*3]));
      v_count++;
      if (post_order) {finalized_vertices[have_finalized] = v_idx; have_finalized++;}
      element_counter++;
      if (element_counter == element_number)
      {
        read_buffer();
      }
      else
      {
        element_descriptor = element_descriptor >> 1;
      }
      return SM_VERTEX;
    }
    else // next element is a triangle
    {
      if (endian_swap) VecCopy3iv_swap_endian(t_idx, (int*)(&element_buffer[element_counter*3]));
      else VecCopy3iv(t_idx, (int*)(&element_buffer[element_counter*3]));
      f_count++;
      for (int i = 0; i < 3; i++)
      {
        if (t_idx[i] < 0)
        {
          t_idx[i] = v_count+t_idx[i];
          t_final[i] = true;
          finalized_vertices[have_finalized] = t_idx[i];
          have_finalized++;
        }
        else
        {
          t_idx[i] = t_idx[i]-1;
          t_final[i] = false;
        }
      }
      element_counter++;
      if (element_counter == element_number)
      {
        read_buffer();
      }
      else
      {
        element_descriptor = element_descriptor >> 1;
      }
      return SM_TRIANGLE;
    }
  }

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

SMevent SMreader_smb::read_event()
{
  if (have_finalized)
  {
    final_idx = finalized_vertices[next_finalized];
    have_finalized--; next_finalized++;
    return SM_FINALIZED;
  }
  else
  {
    return read_element();
  }
}

static int swap_endian_int(int input)
{
  int output;
  ((char*)&output)[0] = ((char*)&input)[3];
  ((char*)&output)[1] = ((char*)&input)[2];
  ((char*)&output)[2] = ((char*)&input)[1];
  ((char*)&output)[3] = ((char*)&input)[0];
  return output;
}

static unsigned int swap_endian_uint(unsigned int input)
{
  int output;
  ((char*)&output)[0] = ((char*)&input)[3];
  ((char*)&output)[1] = ((char*)&input)[2];
  ((char*)&output)[2] = ((char*)&input)[1];
  ((char*)&output)[3] = ((char*)&input)[0];
  return output;
}

#define SM_LITTLE_ENDIAN 0
#define SM_BIG_ENDIAN 1

void SMreader_smb::read_header()
{
  int input;
  // read endianness
#if (defined(i386) || defined(WIN32))   // if little endian machine
  if (fgetc(file) == SM_LITTLE_ENDIAN) endian_swap = false;
  else endian_swap = true;
#else                                   // else big endian machine
  if (fgetc(file) == SM_BIG_ENDIAN) endian_swap = false;
  else endian_swap = true;
#endif
  // read compression flags (not used yet)
  fgetc(file);
  fgetc(file);
  // read comments
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) ncomments = swap_endian_int(input);
  else ncomments = input;
  if (ncomments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      fread(&input, sizeof(int), 1, file);
      if (endian_swap) input = swap_endian_int(input);
      comments[i] = (char*)malloc(sizeof(char)*input);
      fread(comments[i], sizeof(char), input, file);
    }
  }
  // read nverts
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) input = swap_endian_int(input);
  if (input != -1) nverts = input;
  // read nfaces
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) input = swap_endian_int(input);
  if (input != -1) nfaces = input;
  // read bounding box
  if (getc(file))
  {
    if (bb_min_f) delete [] bb_min_f;
    if (bb_max_f) delete [] bb_max_f;
    bb_min_f = new float[3];
    bb_max_f = new float[3];
    if (endian_swap)
    {
      float temp[3];
      fread(temp, sizeof(float), 3, file);
      VecCopy3fv_swap_endian(bb_min_f, temp);
      fread(temp, sizeof(float), 3, file);
      VecCopy3fv_swap_endian(bb_max_f, temp);
    }
    else
    {
      fread(bb_min_f, sizeof(float), 3, file);
      fread(bb_max_f, sizeof(float), 3, file);
    }
  }
}

void SMreader_smb::read_buffer()
{
  fread(&element_descriptor, sizeof(int), 1, file);
  if (endian_swap) element_descriptor = swap_endian_uint(element_descriptor);
  element_number = fread(element_buffer, sizeof(int), 32*3, file) / 3;
  element_counter = 0;
}

SMreader_smb::SMreader_smb()
{
  // init of SMreader_smb
  file = 0;
  have_finalized = 0; next_finalized = 0;

  element_buffer = (int*)malloc(sizeof(int)*3*32);
  element_number = 0;
  element_counter = 0;
}

SMreader_smb::~SMreader_smb()
{
  // clean-up for SMwriter_smb
  free (element_buffer);
}
