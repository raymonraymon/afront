/*
===============================================================================

  FILE:  dynamicqueue.cpp

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
#include "dynamicqueue.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

DynamicQueue::DynamicQueue()
{
  current_capacity = 1024;
  current_capacity_mask = current_capacity-1;
  data = (void**) malloc(sizeof(void*)*current_capacity);
  current_size = 0;
  current_begin = 0;
  current_end = 0;
  number_elements = 0;
}

DynamicQueue::~DynamicQueue()
{
  free(data);
}

int DynamicQueue::size() const
{
  return current_size;
}

int DynamicQueue::elements() const
{
  return number_elements;
}

void* DynamicQueue::getFirstElement() const
{
  return data[current_begin];
}

void* DynamicQueue::getAndRemoveFirstElement()
{
  void* d = data[current_begin];
  removeFirstElement();
  return d;
}

void DynamicQueue::addElement(void* d)
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
      if (data[i]) ((int*)(data[i]))[0] = i; // update indices
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
  number_elements++;
}

void DynamicQueue::removeElement(const void* d)
{
  int i = ((const int*)d)[0];
  if (i == current_begin)
  {
    removeFirstElement();
  }
  else
  {
    data[i] = 0;
    number_elements--;
  }
}

void DynamicQueue::removeFirstElement()
{
  current_begin = (current_begin + 1) & current_capacity_mask;
  current_size--;
  number_elements--;
  while (current_size && (((int*)data)[current_begin] == 0))
  {
    current_begin = (current_begin + 1) & current_capacity_mask;
    current_size--;
  }
}

void* DynamicQueue::iterateGetFirst() const
{
  return data[current_begin];
}

void* DynamicQueue::iterateGetNext(const void* d) const
{
  int i = ((const int*)d)[0];
  while (true)
  {
    i = (i + 1) & current_capacity_mask;
    if (data[i]) return data[i];
  };
}
