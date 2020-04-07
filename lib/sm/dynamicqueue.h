/*
===============================================================================

  FILE:  dynamicqueue.h

  CONTENTS:

    the dynamicqueue class allows adding elements at the end, getting elements
    from the front, but also deleting arbitrary elements in constant time. these
    elements are not really deleted, but marked for deletion. they are actually
    deleted once the element before them was deleted. hence the circular buffer
    used for keeping the elements may increase drastically above the number of
    actually kept elements if the front elements are never removed, but many
    elements are added and only elements in the middle are deleted.

    the elements are expected to provide a field in which their index within the
    dynamic queue can be stored. therefore elements are only stored by references
    (e.g. not by value) in the dynamicqueue. therefore the queue does exclusively
    operate on void* pointers elements.

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2005  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    10 January 2004 -- created after watching 'american sweethearts' with shengi

===============================================================================
*/
#ifndef DYNAMIC_QUEUE_H
#define DYNAMIC_QUEUE_H

class DynamicQueue
{
public:
  DynamicQueue();
  ~DynamicQueue();

  int size() const;
  int elements() const;
  void* getFirstElement() const;
  void* getAndRemoveFirstElement();
  void addElement(void* d);
  void removeElement(const void* d);
  void removeFirstElement();

  void* iterateGetFirst() const;
  void* iterateGetNext(const void* d) const;

private:  
  void** data;
  int current_capacity;
  int current_capacity_mask;
  int current_begin;
  int current_end;
  int current_size;
  int number_elements;
};

#endif
