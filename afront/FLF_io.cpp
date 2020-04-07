
/*

Copyright 2007 University of Utah


This file is part of Afront.

Afront is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Afront is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA

*/


/**
 *
 * file:	FLFio.cpp
 * description:	Feature Loop File converter provides the file
 *		write and read functions that will convert our
 *		extracted features into a format that is necessary
 *		for John and Carlos' triangulation method.
 *
 * author:	Joel Daniels and Linh Ha
 * date:	May 2006
 *
 */

// -- INCLUDES -- //
#include <string.h>
#include "FLF_io.h"

/**
 * write the feature loop file (FLF).  This will output first the list of 
 * points as a set of 'x y z T/F', where the 'T/F' is the true false information
 * indicating if a point is a corner.   This information is then followed by
 * a list of indices for each of the feature loops where each index is coupled
 * with a vector '# x y z'. the indices will start at '0' to indicate the first
 * point in the point list.
 */
void 
FLF_write( const char *filename, 
	   const std::vector< InfoPoint > &ptList,
	   const std::vector< std::vector< InfoIndex > > &featList )
{
  // 0. make sure that we can open the file
  FILE *filePtr = fopen(filename,"w");
  if (!filePtr)
    {
      fprintf(stderr,"[ERROR] FLF_write() couldn't open file '%s'\n",filename);
      return;
    }

  // 1. output some header information to the file
  fprintf(filePtr,"FLF file\n");
  fprintf(filePtr,"# feature loop file automated output\n");
  fprintf(filePtr,"# authors: Joel Daniels and Linh Ha\n");
  fprintf(filePtr,"# date: May 2006\n");
  fprintf(filePtr,"#\n");

  // 2. output the point list
  fprintf(filePtr,"# -- POINT LIST --\n");
  fprintf(filePtr,"%d # number of feature points\n", ptList.size());
  for(unsigned int i=0; i<ptList.size(); i++)
    {
      fprintf(filePtr,"%g %g %g %d # InfoPoint[%d]\n",ptList[i].x_,ptList[i].y_,ptList[i].z_,ptList[i].isCorner_,i);
    }
  fprintf(filePtr,"#\n");

  // 3. output the feature list
  fprintf(filePtr,"# -- FEATURE LIST --\n");
  fprintf(filePtr,"%d # number of features\n",featList.size());
  for(unsigned int i=0; i<featList.size(); i++)
    {
      // output the i'th feature loop
      fprintf(filePtr,"# -- Feature %d --\n",i);
      fprintf(filePtr,"%d # number of feature points in loop\n",featList[i].size());
      for(unsigned int j=0; j<featList[i].size(); j++)
	{
	  fprintf(filePtr,"%d %g %g %g # featurePoint[%d]\n",featList[i][j].index_,featList[i][j].x_,featList[i][j].y_,featList[i][j].z_, i);
	}
    }

  // 4. close the file
  fclose( filePtr );
}

/**
 * read in the FLF file.  This reverses the write function to build
 * the a point list and the feature loop list that indexes the point
 * list.
 */
void 
FLF_read( const char *filename, 
	  std::vector< InfoPoint > &ptList, 
	  std::vector< std::vector< InfoIndex > > &featList )
{
  ptList.clear();
  featList.clear();

  // 0. ensure that we can read in the given file
  std::ifstream in;
  in.open(filename);
  if (!in)
    {
      fprintf(stderr,"[ERROR] FLF_read() couldn't open file: '%s'\n",filename);
      return;
    }

  // 1. make sure that it is really an FLF file
  char buffer[256];
  in.getline( buffer,256 );
  if ( strcmp( buffer,"FLF file" ) != 0 && strcmp( buffer,"FLF file\r" ) != 0 )
    {
      fprintf(stderr,"[ERROR] FLF_read() file missing FLF tag: '%s'\n",filename);
      return;
    }
  fprintf(stderr," *- reading file: '%s'...",filename);

  // 2. read in the header information
  in.getline( buffer,256 );	// # feature loop file automated output
  in.getline( buffer,256 );	// # authors: ...
  in.getline( buffer,256 );	// # date: ...
  in.getline( buffer,256 );	// #

  // 3. read in the point list
  in.getline( buffer,256 );	// # -- POINT LIST --
  unsigned int size;
  in >> size;   
  fprintf(stderr,"reading %d points.",size);
  in.getline( buffer,256 );	// # number of points
  for(unsigned int i=0; i<size; i++)
    {
      // read in a new point information for the point list
      InfoPoint newPoint;
      in >> newPoint.x_ >> newPoint.y_ >> newPoint.z_ >> newPoint.isCorner_;
      ptList.push_back(newPoint);
      in.getline( buffer,256 );	// # InfoPoint[#]
    }
  in.getline( buffer,256 );	// #

  // 4. read in feature list
  in.getline( buffer,256 );	// # -- FEATURE LIST --
  in >> size;
  fprintf(stderr,"reading %d features.",size);
  in.getline( buffer,256 );	// # number of features
  for(unsigned int i=0; i<size; i++)
    {
      in.getline( buffer,256 );	// # -- Feature # --
      unsigned int featureSize;
      in >> featureSize;
      fprintf(stderr,"%d.",featureSize);
      in.getline( buffer,256 );	// # number of feature points in loop
      std::vector< InfoIndex > feature;
      for(unsigned int j=0; j<featureSize; j++)
	{
	  // read in the information for the index and normal
	  InfoIndex newInfo;
	  in >> newInfo.index_ >> newInfo.x_ >> newInfo.y_ >> newInfo.z_;
	  feature.push_back(newInfo);
	  in.getline( buffer,256 );	// # featurePoint[#]
	}
      featList.push_back(feature);
    }

  // 5. close the file
  in.close();
  fprintf(stderr,"done -*\n");
}
