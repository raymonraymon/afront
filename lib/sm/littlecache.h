/*
===============================================================================

  FILE:  littlecache.h

  CONTENTS:


  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2003  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    07 October 2003 -- initial version created the day of the California recall

===============================================================================
*/
#ifndef LITTLE_CACHE_H
#define LITTLE_CACHE_H

class LittleCache
{
public:
  LittleCache(int size = 3);
  ~LittleCache();

  inline void put(void* d0, void* d1, void* d2);
  inline void put(int i0, int i1, int i2, void* d0, void* d1, void* d2);
  inline void put(void* d0, void* d1, void* d2, void* d3);
  inline void* get(int i);
  inline int pos(int i);
  inline int pos(void* d);

//private:
  int size;
  int index[4];
  void* data[4];
};

inline LittleCache::LittleCache(int size)
{
  this->size = size;
  index[0] = -1;
  index[1] = -1;
  index[2] = -1;
  index[3] = -1;
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  data[3] = 0;
}

inline LittleCache::~LittleCache()
{
}

inline void LittleCache::put(void* d0, void* d1, void* d2)
{
  data[0] = d0;
  data[1] = d1;
  data[2] = d2;
}

inline void LittleCache::put(int i0, int i1, int i2, void* d0, void* d1, void* d2)
{
  index[0] = i0;
  index[1] = i1;
  index[2] = i2;
  data[0] = d0;
  data[1] = d1;
  data[2] = d2;
}

inline void LittleCache::put(void* d0, void* d1, void* d2, void* d3)
{
  data[0] = d0;
  data[1] = d1;
  data[2] = d2;
  data[3] = d3;
}

inline void* LittleCache::get(int i)
{
  return data[i];
}

inline int LittleCache::pos(int i)
{
  for (int j = 0; j < size; j++)
  {
    if (i == index[j])
    {
      return j;
    }
  }
  return -1;
}

inline int LittleCache::pos(void* d)
{
  for (int j = 0; j < size; j++)
  {
    if (d == data[j])
    {
      return j;
    }
  }
  return -1;
}

#endif
