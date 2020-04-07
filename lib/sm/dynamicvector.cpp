/*
===============================================================================

  FILE:  dynamicvector.cpp

  CONTENTS:

    the dynamicvector class allows adding and deleting of elements in constant
    time while keeping the n element indexed though numbers 0 to n-1. as both,
    the absolute and the relative index of an element may change over time, the elements are expected to
    provide a field in which their absolute index can be updated. therefore elements are
    only stored by references (e.g. not by value) in the dynamicvector. therefore
    the dynamicvector exclusively operates on void* pointers. 

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2003  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    03 October 2003 -- initial version created the day Peter had a sore throat

===============================================================================
*/
#include "dynamicvector.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

DynamicVector::DynamicVector()
{
  current_capacity = 1024;
  current_capacity_mask = current_capacity-1;
  data = (void**) malloc(sizeof(void*)*current_capacity);
  current_size = 0;
  current_begin = 0;
  current_end = 0;
}

DynamicVector::~DynamicVector()
{
  free(data);
}

int DynamicVector::size() const
{
  return current_size;
}

void* DynamicVector::getFirstElement() const
{
  return data[current_begin];
}

void* DynamicVector::getAndRemoveFirstElement()
{
  if (current_size)
  {
    void* d = data[current_begin];
    current_begin = (current_begin + 1) & current_capacity_mask;
    current_size--;
    return d;
  }
  else
  {
    return 0;
  }
}

void* DynamicVector::getElementWithRelativeIndex(int ri) const
{
  return data[(ri + current_begin) & current_capacity_mask];
}

int DynamicVector::getRelativeIndex(void* d) const
{
  int ai = ((int*)d)[0];
  if (current_begin <= ai)
  {
    return ai - current_begin;
  }
  else
  {
    return ai + current_capacity - current_begin;
  }
}

void DynamicVector::addElement(void* d)
{
  if (current_size == current_capacity)
  {
    void** temp = (void**) malloc(sizeof(void*)*current_capacity*2);
    if (current_begin < current_end)
    {
      memcpy(temp, &(data[current_begin]), sizeof(void*)*current_size);
    }
    else
    {
      int rest = current_size-current_begin;
      memcpy(temp, &(data[current_begin]), sizeof(void*)*rest);
      memcpy(&(temp[rest]), data, sizeof(void*)*(current_size-rest));
    }
    free(data);
    data = temp;
    for (int i = 0; i < current_size; i++)
    {
      ((int*)(data[i]))[0] = i;  // update absolute indices
    }
    current_begin = 0;
    current_end = current_size;
    current_capacity = current_capacity * 2;
    current_capacity_mask = current_capacity - 1;
  }
  data[current_end] = d;
  ((int*)d)[0] = current_end; // set absolute index
  current_end = (current_end + 1) & current_capacity_mask;
  current_size++;
}

void DynamicVector::removeElement(void* d)
{
  int ai = ((int*)d)[0];
  if (ai != current_begin) // *not* removing the first element (e.g. the one at current_begin)
  {
    data[ai] = data[current_begin];                          // move first element to where the deleted element used to be
    ((int*)(data[ai]))[0] = ai;                              // update absolute index of first element
  }
  current_begin = (current_begin + 1) & current_capacity_mask;
  current_size--;
}

void DynamicVector::removeFirstElement()
{
  if (current_size)
  {
    current_begin = (current_begin + 1) & current_capacity_mask;
    current_size--;
  }
}
