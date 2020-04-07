/*
===============================================================================

  FILE:  smreader_sma.cpp
  
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
#include "smreader_sma.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "vec3fv.h"
#include "vec3iv.h"

#define MAX_LINE_SIZE 512

bool SMreader_sma::open(FILE* file)
{
  if (file == 0)
  {
    return false;
  }

#ifdef _WIN32
  if (file == stdin)
  {
    if(_setmode( _fileno( stdin ), _O_TEXT ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdin to text (translated) mode\n");
    }
  }
#endif

  this->file = file;

  skipped_lines = 0;
  line = (char*)malloc(sizeof(char)*MAX_LINE_SIZE);
  if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
  {
    free(line);
    line = 0;
    return false;
  }

  // look for header information
  char dummy[MAX_LINE_SIZE];

  while (line && line[0] != 'v' && line[0] != 'f' && line[0] != 'x')
  {
    if (strstr(line, "nverts"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &nverts);
    }
    else if (strstr(line, "nfaces"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &nfaces);
    }
    else if (strstr(line, "bb_min"))
    {
      if (bb_min_f == 0) bb_min_f = new float[3];
      sscanf(&(line[1]), "%s %f %f %f", dummy, &(bb_min_f[0]), &(bb_min_f[1]), &(bb_min_f[2]));
    }
    else if (strstr(line, "bb_max"))
    {
      if (bb_max_f == 0) bb_max_f = new float[3];
      sscanf(&(line[1]), "%s %f %f %f", dummy, &(bb_max_f[0]), &(bb_max_f[1]), &(bb_max_f[2]));
    }
    else if (line[0] == '#')
    {
      if (sscanf(&(line[1]), "%s", dummy))
      {
        if (comments == 0)
        {
          ncomments = 0;
          comments = (char**)malloc(sizeof(char*)*10);
          comments[9] = (char*)-1;
        }
        else if (comments[ncomments] == (char*)-1)
        {
          comments = (char**)realloc(comments,sizeof(char*)*ncomments*2);
          comments[ncomments*2-1] = (char*)-1;
        }
        comments[ncomments] = strdup(dummy);
        ncomments++;
      }
    }
    else
    {
      int j = (int)strlen(line)-1;
      if (j > 0)
      {
        for (int i = 0; i < j; i++)
        {
          if (!isprint(line[i]))
          {
            fprintf(stderr,"FATAL ERROR: input file seems to be binary\n");
            break;
          }
        }
        line[(10<j?10:j)] = '\0';
        if (skipped_lines < 5)
        {
          fprintf(stderr,"WARNING: skipping line '%s...'.\n",line);
        }
        else if (skipped_lines == 5)
        {
          fprintf(stderr,"WARNING: skipping more than 5 lines. being quiet.\n");
        }
        skipped_lines++;
      }
    }
    if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
    {
      free(line);
      line = 0;
      return false;
    }
  }
  if (line[0] == 'v')
  {
    post_order = false;
  }
  else if (line[0] == 'f')
  {
    post_order = true;
  }
  else
  {
    fprintf(stderr,"WARNING: cannot determine order. assuming pre-order SM.\n");
  }

  v_count = 0;
  f_count = 0;

  return true;
}

void SMreader_sma::close(bool close_file)
{
  // close of SMreader_sma
  if (skipped_lines) fprintf(stderr,"WARNING: skipped %d lines.\n",skipped_lines);
  skipped_lines = 0;
  if (close_file && file /*&& file != stdin*/) fclose(file);
  file = 0;
  if (line)
  {
    free(line);
    line = 0;
  }
  have_finalized = 0; next_finalized = 0;

  // close of SMreader interface
  v_count = -1;
  f_count = -1;
}

SMevent SMreader_sma::read_element()
{
  have_finalized = next_finalized = 0;
  while (line)
  {
    if ((line[0] == 'v') && (line[1] == ' '))
    {
      v_idx = v_count;
      sscanf(&(line[1]), "%f %f %f", &(v_pos_f[0]), &(v_pos_f[1]), &(v_pos_f[2]));
      if (post_order) {finalized_vertices[have_finalized] = v_idx; have_finalized++;}
      v_count++;
      if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SM_VERTEX;
    }
    else if ((line[0] == 'f') && (line[1] == ' '))
    {
      sscanf(&(line[1]), "%d %d %d", &(t_idx[0]), &(t_idx[1]), &(t_idx[2]));
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
      if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SM_TRIANGLE;
    }
    else if ((line[0] == 'x') && (line[1] == ' '))
    {
      sscanf(&(line[1]), "%d", &(final_idx));
      if (final_idx < 0)
      {
        final_idx = v_count+final_idx;
      }
      else
      {
        final_idx = final_idx - 1;
      }
      if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SM_FINALIZED;
    }
    else if (line[0] == '#')
    {
      // comments in the body are silently ignored
      if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
      {
        free(line);
        line = 0;
      }
    }
    else
    {
      int j = (int)strlen(line)-1;
      if (j)
      {
        for (int i = 0; i < j; i++)
        {
          if (!isprint(line[i]))
          {
            fprintf(stderr,"WARNING: input file seems to be binary\n");
            break;
          }
        }
        line[(10<j?10:j)] = '\0';
        if (skipped_lines < 5)
        {
          fprintf(stderr,"WARNING: skipping line '%s...'.\n",line);
        }
        else if (skipped_lines == 5)
        {
          fprintf(stderr,"WARNING: skipping more than 5 lines. being quiet.\n");
        }
        skipped_lines++;
      }
      if (fgets(line, sizeof(char) * MAX_LINE_SIZE, file) == 0)
      {
        free(line);
        line = 0;
      }
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

SMevent SMreader_sma::read_event()
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

SMreader_sma::SMreader_sma()
{
  // init of SMreader_sma
  file = 0;
  line = 0;
  have_finalized = 0; next_finalized = 0;
}

SMreader_sma::~SMreader_sma()
{
  // clean-up for SMreader_sma
}
