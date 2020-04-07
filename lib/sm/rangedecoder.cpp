/*
===============================================================================

  FILE:  rangedecoder.cpp
  
  CONTENTS:
      
    see header file

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003 martin isenburg (isenburg@cs.unc.edu)
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see header file
  
===============================================================================
*/
#include "rangedecoder.h"

inline void RangeDecoder::normalize()
{
  while (range <= BOTTOM_VALUE)
  {
    low = (low<<8) | ((buffer<<EXTRA_BITS)&0xff);
    buffer = inbyte();
    low |= buffer >> (8-EXTRA_BITS);
    range <<= 8;
  }
}

/* Decode with modelling                                     */
unsigned int RangeDecoder::decode(RangeModel* rm)
{
  unsigned int sym;
  unsigned int ltfreq;
  unsigned int syfreq;
  unsigned int tmp;
  unsigned int lg_totf = rm->lg_totf;

  normalize();
  help = this->range>>lg_totf;
  ltfreq = low/help;
#ifdef EXTRAFAST
  ltfreq = ltfreq;
#else
  ltfreq = ((ltfreq>>lg_totf) ? (1<<lg_totf)-1 : ltfreq);
#endif

  sym = rm->getsym(ltfreq);
  rm->getfreq(sym,&syfreq,&ltfreq);

  tmp = help * ltfreq;
  low -= tmp;
#ifdef EXTRAFAST
  this->range = help * syfreq;
#else
  if ((ltfreq + syfreq) < (unsigned int)(1<<lg_totf))
  {
    this->range = help * syfreq;
  }
  else
  {
    this->range -= tmp;
  }
#endif

  rm->update(sym);

  return sym;
}

/* Decode a range without modelling                          */
unsigned int RangeDecoder::decode(unsigned int range)
{
  unsigned int tmp;
  unsigned int tmp1;

  if (range > 4194303) // 22 bits
  {
    tmp = decodeShort();
    range = range >> 16;
    range++;
    tmp1 = decode(range) << 16;
    return (tmp1|tmp);
  }
  
  normalize();
  help = this->range/range;
  tmp = low/help;
#ifdef EXTRAFAST
  tmp = tmp;
#else
  tmp = (tmp>=range ? range-1 : tmp);
#endif

  tmp1 = (help * tmp);
  low -= tmp1;
#ifdef EXTRAFAST
  this->range = help;
#else
  if (tmp+1 < range)
  {
    this->range = help;
  }
  else
  {
    this->range -= tmp1;
  }
#endif

  return tmp;
}

/* Decode a range without modelling                          */
unsigned INT64 RangeDecoder::decode64(unsigned INT64 range)
{
  if (range > 4294967295u) // 32 bits
  {
    unsigned INT64 tmp;
    unsigned INT64 tmp1;
    tmp = decodeInt();
    range = range >> 32;
    range++;
    tmp1 = ((INT64)decode((unsigned int)(range))) << 32;
    return (tmp1|tmp);
  }
  else
  {
    return decode((unsigned int)(range));
  }
}

/* Decode a byte without modelling                           */
unsigned char RangeDecoder::decodeByte()
{
  unsigned char tmp = culshift(8);
  update(1, tmp, (unsigned int)1<<8);
  return tmp;
}

/* Decode a short without modelling                          */
unsigned short RangeDecoder::decodeShort()
{
  unsigned short tmp = culshift(16);
  update(1, tmp, (unsigned int)1<<16);
  return tmp;
}

/* Decode an unsigned int without modelling                  */
unsigned int RangeDecoder::decodeInt()
{
  unsigned int lowerInt = decodeShort();
  unsigned int upperInt = decodeShort();
  return upperInt*65536+lowerInt;
}

/* Decode a float without modelling                          */
float RangeDecoder::decodeFloat()
{
  float f;
  *((unsigned int*)(&f)) = decodeInt();
  return f;
}

/* Decode an unsigned int64 without modelling                */
unsigned INT64 RangeDecoder::decodeInt64()
{
  unsigned INT64 lowerInt = decodeInt();
  unsigned INT64 upperInt = decodeInt();
  return upperInt*4294967296ll+lowerInt;
}

/* Decode a double without modelling                         */
double RangeDecoder::decodeDouble()
{
  double d;
  *((unsigned INT64*)(&d)) = decodeInt64();
  return d;
}

/* Finish decoding                                           */
void RangeDecoder::done()
{
  normalize();      /* normalize to use up all bytes */
}

unsigned int RangeDecoder::culshift(unsigned int shift)
{
  unsigned int tmp;
  normalize();
  help = range>>shift;
  tmp = low/help;
#ifdef EXTRAFAST
  return tmp;
#else
  return (tmp>>shift ? ((unsigned int)1<<shift)-1 : tmp);
#endif
}

/* Update decoding state                                     */
/* sy_f is the interval length (frequency of the symbol)     */
/* lt_f is the lower end (frequency sum of < symbols)        */
/* tot_f is the total interval length (total frequency sum)  */
void RangeDecoder::update(unsigned int sy_f, unsigned int lt_f, unsigned int tot_f)
{
  unsigned int tmp;
  tmp = help * lt_f;
  low -= tmp;
#ifdef EXTRAFAST
  this->range = help * sy_f;
#else
  if (lt_f + sy_f < tot_f)
  {
    this->range = help * sy_f;
  }
  else
  {
    this->range -= tmp;
  }
#endif
}

RangeDecoder::RangeDecoder(unsigned char* chars, int number_chars)
{
  this->chars = chars;
  this->number_chars = number_chars;
  current_char = 0;
  fp = 0;

  buffer = inbyte();
  if (buffer != HEADERBYTE)
  {
    fprintf(stderr, "RangeDecoder: wrong HEADERBYTE of %d. is should be %d\n", buffer, HEADERBYTE);
    return;
  }
  buffer = inbyte();
  low = buffer >> (8-EXTRA_BITS);
  range = (unsigned int)1 << EXTRA_BITS;
}

RangeDecoder::RangeDecoder(FILE* fp)
{
  chars = 0;
  number_chars = 0;
  current_char = 0;
  this->fp = fp;

  buffer = inbyte();
  if (buffer != HEADERBYTE)
  {
    fprintf(stderr, "RangeDecoder: wrong HEADERBYTE of %d. is should be %d\n", buffer, HEADERBYTE);
    return;
  }
  buffer = inbyte();
  low = buffer >> (8-EXTRA_BITS);
  range = (unsigned int)1 << EXTRA_BITS;
}

RangeDecoder::~RangeDecoder()
{
  if (fp)
  {
    fclose(fp);
  }
}

inline unsigned int RangeDecoder::inbyte()
{
  int c;
  if (fp)
  {
    c = getc(fp);
  }
  else
  {
    if (current_char < number_chars)
    {
      c = chars[current_char++];
    }
    else
    {
      c = EOF;
    }
  }
  return c;
}
