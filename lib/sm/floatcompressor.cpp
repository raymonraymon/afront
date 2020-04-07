/*
===============================================================================

  FILE:  floatcompressor.cpp
  
  CONTENTS:
  
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    Copyright (C) 2000  Martin Isenburg (isenburg@cs.unc.edu)
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    30 September 2003 -- created initial version after Ajith's good-bye lunch
  
===============================================================================
*/

#include <stdlib.h>
#include "floatcompressor.h"


// 0  1  2  3  4   5   6   7    8    9  10     11    12  13  14  15  16  17  18   19   20   21    22    23    24

const int bitsLowTable[] = {
   0, 0, 1, 2, 3,  4,  5,  6,   7,   8,  9,    10,   11,  1,  2,  3,  4,  5,  6,   7,   8,   9,   10,   11,   12
};

const int bitsMaskTable[] = {
   0, 0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047,  1,  3,  7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4096
};

const int bitsHighTable[] = {
   0, 0, 0, 0, 0,  0,  0,  0,   0,   0,   0,    0,    0, 11, 11, 11, 11, 11, 11,  11,  11,  11,   11,   11,   11
};

FloatCompressor::FloatCompressor()
{
  bits = 0;
  bits_for_mantissa = (int*)malloc(sizeof(int)*256);
}

FloatCompressor::~FloatCompressor()
{
  free(bits_for_mantissa);
}

int num_predictions_none = 0;
int num_predictions_last = 0;
int num_predictions_across = 0;
int num_predictions_within = 0;

void FloatCompressor::SetupCompressor(RangeEncoder* re, bool within)
{
  if (re == 0)
  {
    // none predictions
    ae_sign_none = new RangeEncoder(0);
    ae_exponent_none = new RangeEncoder(0);
    ae_mantissa_none = new RangeEncoder(0);

    // last predictions
    ae_sign_last = new RangeEncoder(0);
    ae_exponent_last = new RangeEncoder(0);
    ae_mantissa_last = new RangeEncoder(0);

    // across predictions
    ae_sign_across = new RangeEncoder(0);
    ae_exponent_across = new RangeEncoder(0);
    ae_mantissa_across = new RangeEncoder(0);

    if (within)
    {
      // use within predictions
      ae_sign_within = new RangeEncoder(0);
      ae_exponent_within = new RangeEncoder(0);
      ae_mantissa_within = new RangeEncoder(0);
    }
    else
    {
      // do not use within predictions
      ae_sign_within = 0;
      ae_exponent_within = 0;
      ae_mantissa_within = 0;
    }
  }
  else
  {
    // none predictions
    ae_sign_none = re;
    ae_exponent_none = re;
    ae_mantissa_none = re;

    // last predictions
    ae_sign_last = re;
    ae_exponent_last = re;
    ae_mantissa_last = re;

    // across predictions
    ae_sign_across = re;
    ae_exponent_across = re;
    ae_mantissa_across = re;

    if (within)
    {
      // use within predictions
      ae_sign_within = re;
      ae_exponent_within = re;
      ae_mantissa_within = re;
    }
    else
    {
      // do not use within predictions
      ae_sign_within = 0;
      ae_exponent_within = 0;
      ae_mantissa_within = 0;
    }
  }
  alloc_range_tables(within);
  reset_precision();
}

void FloatCompressor::SetupDecompressor(RangeDecoder* rd, bool within)
{
  if (rd == 0)
  {
    // none predictions
    ad_sign_none = new RangeDecoder(ae_sign_none->getChars(), ae_sign_none->getNumberChars());
    ad_exponent_none = new RangeDecoder(ae_exponent_none->getChars(), ae_exponent_none->getNumberChars());
    ad_mantissa_none = new RangeDecoder(ae_mantissa_none->getChars(), ae_mantissa_none->getNumberChars());

    // last predictions
    ad_sign_last = new RangeDecoder(ae_sign_last->getChars(), ae_sign_last->getNumberChars());
    ad_exponent_last = new RangeDecoder(ae_exponent_last->getChars(), ae_exponent_last->getNumberChars());
    ad_mantissa_last = new RangeDecoder(ae_mantissa_last->getChars(), ae_mantissa_last->getNumberChars());

    // across predictions
    ad_sign_across = new RangeDecoder(ae_sign_across->getChars(), ae_sign_across->getNumberChars());
    ad_exponent_across = new RangeDecoder(ae_exponent_across->getChars(), ae_exponent_across->getNumberChars());
    ad_mantissa_across = new RangeDecoder(ae_mantissa_across->getChars(), ae_mantissa_across->getNumberChars());

    if (within)
    {
      // use within predictions
      ad_sign_within = new RangeDecoder(ae_sign_within->getChars(), ae_sign_within->getNumberChars());
      ad_exponent_within = new RangeDecoder(ae_exponent_within->getChars(), ae_exponent_within->getNumberChars());
      ad_mantissa_within = new RangeDecoder(ae_mantissa_within->getChars(), ae_mantissa_within->getNumberChars());
    }
    else
    {
      // do not use within predictions
      ad_sign_within = 0;
      ad_exponent_within = 0;
      ad_mantissa_within = 0;
    }
  }
  else
  {
    // none predictions
    ad_sign_none = rd;
    ad_exponent_none = rd;
    ad_mantissa_none = rd;

    // last predictions
    ad_sign_last = rd;
    ad_exponent_last = rd;
    ad_mantissa_last = rd;

    // across predictions
    ad_sign_across = rd;
    ad_exponent_across = rd;
    ad_mantissa_across = rd;

    if (within)
    {
      // use within predictions
      ad_sign_within = rd;
      ad_exponent_within = rd;
      ad_mantissa_within = rd;
    }
    else
    {
      // do not use within predictions
      ad_sign_within = 0;
      ad_exponent_within = 0;
      ad_mantissa_within = 0;
    }
  }
  alloc_range_tables(within);
  reset_precision();
}

void FloatCompressor::alloc_range_tables(bool within)
{
  int i;

  rmSignNone = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmExponentNone = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmMantissaNoneLow = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmMantissaNoneHigh = (RangeModel**)malloc(sizeof(RangeModel*)*256);

  rmSignLast = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmExponentLast = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmMantissaLastLow = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmMantissaLastHigh = (RangeModel**)malloc(sizeof(RangeModel*)*256);

  rmSignAcross = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmExponentAcross = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmMantissaAcrossLow = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmMantissaAcrossHigh = (RangeModel**)malloc(sizeof(RangeModel*)*256);

  if (within)
  {
    rmSignWithin = (RangeModel**)malloc(sizeof(RangeModel*)*256);
    rmExponentWithin = (RangeModel**)malloc(sizeof(RangeModel*)*256);
    rmMantissaWithinLow = (RangeModel**)malloc(sizeof(RangeModel*)*256);
    rmMantissaWithinHigh = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  }
  else
  {
    rmSignWithin = 0;
    rmExponentWithin = 0;
    rmMantissaWithinLow = 0;
    rmMantissaWithinHigh = 0;
  }

  rmZeroOrExtremeValue = (RangeModel**)malloc(sizeof(RangeModel*)*256);
  rmPositiveOrNegativeValue = (RangeModel**)malloc(sizeof(RangeModel*)*256);

  for (i = 0; i < 256; i++)
  {
    rmSignNone[i] = 0;
    rmExponentNone[i] = 0;
    rmMantissaNoneLow[i] = 0;
    rmMantissaNoneHigh[i] = 0;

    rmSignLast[i] = 0;
    rmExponentLast[i] = 0;
    rmMantissaLastLow[i] = 0;
    rmMantissaLastHigh[i] = 0;

    rmSignAcross[i] = 0;
    rmExponentAcross[i] = 0;
    rmMantissaAcrossLow[i] = 0;
    rmMantissaAcrossHigh[i] = 0;

    if (within)
    {
      rmSignWithin[i] = 0;
      rmExponentWithin[i] = 0;
      rmMantissaWithinLow[i] = 0;
      rmMantissaWithinHigh[i] = 0;
    }

    rmZeroOrExtremeValue[i] = 0;
    rmPositiveOrNegativeValue[i] = 0;
  }
}

void FloatCompressor::dealloc_range_tables(bool within)
{
  int i;

  for (i = 0; i < 256; i++)
  {
    if (rmSignNone[i]) delete rmSignNone[i];
    if (rmExponentNone[i]) delete rmExponentNone[i];
    if (rmMantissaNoneLow[i]) delete rmMantissaNoneLow[i];
    if (rmMantissaNoneHigh[i]) delete rmMantissaNoneHigh[i];

    if (rmSignLast[i]) delete rmSignLast[i];
    if (rmExponentLast[i]) delete rmExponentLast[i];
    if (rmMantissaLastLow[i]) delete rmMantissaLastLow[i];
    if (rmMantissaLastHigh[i]) delete rmMantissaLastHigh[i];

    if (rmSignAcross[i]) delete rmSignAcross[i];
    if (rmExponentAcross[i]) delete rmExponentAcross[i];
    if (rmMantissaAcrossLow[i]) delete rmMantissaAcrossLow[i];
    if (rmMantissaAcrossHigh[i]) delete rmMantissaAcrossHigh[i];

    if (within)
    {
      if (rmSignWithin[i]) delete rmSignWithin[i];
      if (rmExponentWithin[i]) delete rmExponentWithin[i];
      if (rmMantissaWithinLow[i]) delete rmMantissaWithinLow[i];
      if (rmMantissaWithinHigh[i]) delete rmMantissaWithinHigh[i];
    }

    if (rmZeroOrExtremeValue[i]) delete rmZeroOrExtremeValue[i];
    if (rmPositiveOrNegativeValue[i]) delete rmPositiveOrNegativeValue[i];
  }

  free(rmSignNone);
  free(rmExponentNone);
  free(rmMantissaNoneLow);
  free(rmMantissaNoneHigh);

  free(rmSignLast);
  free(rmExponentLast);
  free(rmMantissaLastLow);
  free(rmMantissaLastHigh);

  free(rmSignAcross);
  free(rmExponentAcross);
  free(rmMantissaAcrossLow);
  free(rmMantissaAcrossHigh);

  if (within)
  {
    free(rmSignWithin);
    free(rmExponentWithin);
    free(rmMantissaWithinLow);
    free(rmMantissaWithinHigh);
  }

  free(rmZeroOrExtremeValue);
  free(rmPositiveOrNegativeValue);
}

void FloatCompressor::FinishCompressor(bool within)
{
  if (ae_sign_none != ae_sign_last)
  {
    fprintf(stderr,"none %d last %d across %d within %d\n", num_predictions_none, num_predictions_last, num_predictions_across, num_predictions_within);
    ae_sign_none->done();
    ae_sign_last->done();
    ae_sign_across->done();
    fprintf(stderr,"TOTAL: sign bytes %d bpf %5.2f bpv %5.2f\n",ae_sign_none->getNumberChars()+ae_sign_last->getNumberChars()+ae_sign_across->getNumberChars(),
      (float)(ae_sign_none->getNumberBits()+ae_sign_last->getNumberBits()+ae_sign_across->getNumberBits())/(float)(num_predictions_none+num_predictions_last+num_predictions_across),
      (float)(ae_sign_none->getNumberBits()+ae_sign_last->getNumberBits()+ae_sign_across->getNumberBits())/(float)(num_predictions_none+num_predictions_last+num_predictions_across)*3.0f);
    ae_exponent_none->done();
    ae_exponent_last->done();
    ae_exponent_across->done();
    fprintf(stderr,"TOTAL: expo bytes %d bpf %5.2f bpv %5.2f\n",ae_exponent_none->getNumberChars()+ae_exponent_last->getNumberChars()+ae_exponent_across->getNumberChars(),
      (float)(ae_exponent_none->getNumberBits()+ae_exponent_last->getNumberBits()+ae_exponent_across->getNumberBits())/(float)(num_predictions_none+num_predictions_last+num_predictions_across),
      (float)(ae_exponent_none->getNumberBits()+ae_exponent_last->getNumberBits()+ae_exponent_across->getNumberBits())/(float)(num_predictions_none+num_predictions_last+num_predictions_across)*3.0f);
    ae_mantissa_none->done();
    ae_mantissa_last->done();
    ae_mantissa_across->done();
    fprintf(stderr,"TOTAL: mant bytes %d bpf %5.2f bpv %5.2f\n",ae_mantissa_none->getNumberChars()+ae_mantissa_last->getNumberChars()+ae_mantissa_across->getNumberChars(),
      (float)(ae_mantissa_none->getNumberBits()+ae_mantissa_last->getNumberBits()+ae_mantissa_across->getNumberBits())/(float)(num_predictions_none+num_predictions_last+num_predictions_across),
      (float)(ae_mantissa_none->getNumberBits()+ae_mantissa_last->getNumberBits()+ae_mantissa_across->getNumberBits())/(float)(num_predictions_none+num_predictions_last+num_predictions_across)*3.0f);

    if (ae_sign_none == ae_exponent_none)
    {
      fprintf(stderr,"none: total bytes %d\n",ae_sign_none->getNumberChars());
    }
    else
    {
      fprintf(stderr,"none: sign bytes %d bpf %5.2f bpv %5.2f\n",ae_sign_none->getNumberChars(),(float)ae_sign_none->getNumberBits()/(float)num_predictions_none,(float)ae_sign_none->getNumberBits()/(float)num_predictions_none*3.0f);
      fprintf(stderr,"none: expo bytes %d bpf %5.2f bpv %5.2f\n",ae_exponent_none->getNumberChars(),(float)ae_exponent_none->getNumberBits()/(float)num_predictions_none,(float)ae_exponent_none->getNumberBits()/(float)num_predictions_none*3.0f);
      fprintf(stderr,"none: mant bytes %d bpf %5.2f bpv %5.2f\n",ae_mantissa_none->getNumberChars(),(float)ae_mantissa_none->getNumberBits()/(float)num_predictions_none,(float)ae_mantissa_none->getNumberBits()/(float)num_predictions_none*3.0f);
    }

    if (ae_sign_last == ae_exponent_last)
    {
      fprintf(stderr,"Last: total bytes %d\n",ae_sign_last->getNumberChars());
    }
    else
    {
      fprintf(stderr,"Last: sign bytes %d bpf %5.2f bpv %5.2f\n",ae_sign_last->getNumberChars(),(float)ae_sign_last->getNumberBits()/(float)num_predictions_last,(float)ae_sign_last->getNumberBits()/(float)num_predictions_last*3.0f);
      fprintf(stderr,"Last: expo bytes %d bpf %5.2f bpv %5.2f\n",ae_exponent_last->getNumberChars(),(float)ae_exponent_last->getNumberBits()/(float)num_predictions_last,(float)ae_exponent_last->getNumberBits()/(float)num_predictions_last*3.0f);
      fprintf(stderr,"Last: mant bytes %d bpf %5.2f bpv %5.2f\n",ae_mantissa_last->getNumberChars(),(float)ae_mantissa_last->getNumberBits()/(float)num_predictions_last,(float)ae_mantissa_last->getNumberBits()/(float)num_predictions_last*3.0f);
    }

    if (ae_sign_across == ae_exponent_across)
    {
      fprintf(stderr,"Across: total bytes %d\n",ae_sign_across->getNumberChars());
    }
    else
    {
      fprintf(stderr,"Across: sign bytes %d bpf %5.2f bpv %5.2f\n",ae_sign_across->getNumberChars(),(float)ae_sign_across->getNumberBits()/(float)num_predictions_across,(float)ae_sign_across->getNumberBits()/(float)num_predictions_across*3.0f);
      fprintf(stderr,"Across: expo bytes %d bpf %5.2f bpv %5.2f\n",ae_exponent_across->getNumberChars(),(float)ae_exponent_across->getNumberBits()/(float)num_predictions_across,(float)ae_exponent_across->getNumberBits()/(float)num_predictions_across*3.0f);
      fprintf(stderr,"Across: mant bytes %d bpf %5.2f bpv %5.2f\n",ae_mantissa_across->getNumberChars(),(float)ae_mantissa_across->getNumberBits()/(float)num_predictions_across,(float)ae_mantissa_across->getNumberBits()/(float)num_predictions_across*3.0f);
    }
  }
  dealloc_range_tables(within);
}

void FloatCompressor::FinishDecompressor(bool within)
{
  if (ad_sign_none != ad_sign_last)
  {
    ad_sign_none->done();
    ad_sign_last->done();
    ad_sign_across->done();

    ad_exponent_none->done();
    ad_exponent_last->done();
    ad_exponent_across->done();

    ad_mantissa_none->done();
    ae_mantissa_last->done();
    ae_mantissa_across->done();
  }
  dealloc_range_tables(within);
}

void FloatCompressor::compress_sign(int exponent, int sign, RangeEncoder* re_sign, RangeModel** rmSign)
{
  if (rmSign[exponent] == 0)
  {
    rmSign[exponent] = new RangeModel(2, 0, 1);
  }
  re_sign->encode(rmSign[exponent], sign);
//  fprintf(stderr,"compress_sign exp %d sign %d\n",exponent,sign);
}

int FloatCompressor::decompress_sign(int exponent, RangeDecoder* rd_sign, RangeModel** rmSign)
{
  if (rmSign[exponent] == 0)
  {
    rmSign[exponent] = new RangeModel(2, 0, 0);
  }
  int sign = rd_sign->decode(rmSign[exponent]);
//  fprintf(stderr,"decompress_sign exp %d sign %d\n",exponent,sign);
  return sign;
}

void FloatCompressor::compress_exponent(int exponentPred, int exponentReal, RangeEncoder* re_exponent, RangeModel** rmExponent)
{
  if (rmExponent[exponentPred] == 0)
  {
    rmExponent[exponentPred] = new RangeModel(256,0,1);
  }
  re_exponent->encode(rmExponent[exponentPred],exponentReal);
//  fprintf(stderr,"compress_exponent exp %d real %d\n",exponentPred,exponentReal);
}

int FloatCompressor::decompress_exponent(int exponentPred, RangeDecoder* rd_exponent, RangeModel** rmExponent)
{
  if (rmExponent[exponentPred] == 0)
  {
    rmExponent[exponentPred] = new RangeModel(256,0,0);
  }

  int exponentReal = rd_exponent->decode(rmExponent[exponentPred]);
//  fprintf(stderr,"decompress_exponent exp %d real %d\n",exponentPred,exponentReal);
  return exponentReal;
}

int FloatCompressor::compress_mantissa(int exponent, int mantissaPred, int mantissaReal, RangeEncoder* re_mantissa, RangeModel** rmMantissaHigh, RangeModel** rmMantissaLow)
{
  int c, sign;

  int bits_used = bits_for_mantissa[exponent];

  if (!bits_used) // all mantissa bits are discarded
  {
    return 0;
  }

  int bits_discarded = 23 - bits_used;

  int bits_mask = (1<<23)-(1<<bits_discarded);

  mantissaReal &= bits_mask;

  c = mantissaReal - (mantissaPred&bits_mask);

  if (c < (1-(1<<22)))
  {
    c += (1<<23);
  }
  else if (c > (1<<22))
  {
    c -= (1<<23); 
  }

  // now c should lie in the interval [ - (2^22 - 1)  ...  + (2^22) ]

  if (c == (1<<22))            // the extreme misprediction
  {
    sign = -2;
    c = 0;
  }
  else
  {
    c = c >> bits_discarded;  // discard the zero bits

    if (c > 0)                // positive
    {
      sign = 0;
    }
    else if (c < 0)            // negative
    {
      sign = 1;
      c = -c;
    }
    else                      // zero
    {
      sign = -1;
    }
  }

  int bitsLowUsed = bitsLowTable[bits_used];

  if (bitsLowUsed)
  {
    int bitsHighUsed = bitsHighTable[bits_used];

    if (rmMantissaLow[exponent] == 0)
    {
      rmMantissaLow[exponent] = new RangeModel((1<<bitsLowUsed),0,1);
      if (bitsHighUsed)
      {
        rmMantissaHigh[exponent] = new RangeModel((1<<bitsHighUsed),0,1);
      }
    }

    if (bitsHighUsed)
    {
      int bitsHigh = c >> bitsLowUsed;
      re_mantissa->encode(rmMantissaHigh[exponent], bitsHigh);
      int bitsLow = c & bitsMaskTable[bits_used];
      re_mantissa->encode(rmMantissaLow[exponent], bitsLow);
    }
    else
    {
      re_mantissa->encode(rmMantissaLow[exponent], c);
    }
  }

  if (sign < 0)      // zero or extreme
  {
    if (rmZeroOrExtremeValue[exponent] == 0)
    {
      rmZeroOrExtremeValue[exponent] = new RangeModel(2,0,1);
    }
    if (sign == -1) // zero
    {
      re_mantissa->encode(rmZeroOrExtremeValue[exponent], 0);
    }
    else            // extreme
    {
      re_mantissa->encode(rmZeroOrExtremeValue[exponent], 1);
    }
  }
  else              // positive or negative
  {
    if (rmPositiveOrNegativeValue[exponent] == 0)
    {
      rmPositiveOrNegativeValue[exponent] = new RangeModel(2,0,1);
    }
    if (sign == 0)  // positive
    {
      re_mantissa->encode(rmPositiveOrNegativeValue[exponent], 0);
    }
    else            // negative
    {
      re_mantissa->encode(rmPositiveOrNegativeValue[exponent], 1);
    }
  }
  return mantissaReal;
}

int FloatCompressor::decompress_mantissa(int exponent, int mantissaPred, RangeDecoder* rd_mantissa, RangeModel** rmMantissaLow, RangeModel** rmMantissaHigh)
{
  int c = 0;

  int bits_used = bits_for_mantissa[exponent];

  if (!bits_used) // all mantissa bits are discarded
  {
    return 0;
  }

  int bitsLowUsed = bitsLowTable[bits_used];

  if (bitsLowUsed)
  {
    int bitsHighUsed = bitsHighTable[bits_used];

    if (rmMantissaLow[exponent] == 0)
    {
      rmMantissaLow[exponent] = new RangeModel((1<<bitsLowUsed),0,0);
      if (bitsHighUsed)
      {
        rmMantissaHigh[exponent] = new RangeModel((1<<bitsHighUsed),0,0);
      }
    }

    if (bitsHighUsed)
    {
      c = (rd_mantissa->decode(rmMantissaHigh[exponent]) << bitsLowUsed);
      c = c | rd_mantissa->decode(rmMantissaLow[exponent]);
    }
    else
    {
      c = rd_mantissa->decode(rmMantissaLow[exponent]);
    }
  }

  int bits_discarded = 23 - bits_used;

  if (c != 0) // decode sign
  {
    if (rmPositiveOrNegativeValue[exponent] == 0)
    {
      rmPositiveOrNegativeValue[exponent] = new RangeModel(2,0,0);
    }
    if (rd_mantissa->decode(rmPositiveOrNegativeValue[exponent]))
    {
      c = -(c << bits_discarded);
    }
    else
    {
      c = c << bits_discarded;
    }
  }
  else // decode zero or extreme
  {
    if (rmZeroOrExtremeValue[exponent] == 0)
    {
      rmZeroOrExtremeValue[exponent] = new RangeModel(2,0,0);
    }
    if (rd_mantissa->decode(rmZeroOrExtremeValue[exponent]))
    {
      c = (1<<(22));
    }
  }

  int bits_mask = (1<<23)-(1<<bits_discarded);

  c += (mantissaPred&bits_mask);

  // wrap back into valid range
  if (c < 0)
  {
    c += (1<<23);
  }
  else if (c >= (1<<23))
  {
    c -= (1<<23);
  }

  return c;
}

void FloatCompressor::reset_precision()
{
  range_exponent = -2;
  for (int i = 0; i < 256; i++)
  {
    bits_for_mantissa[i] = 23;
  }
}

void FloatCompressor::update_precision(float number)
{
  if (range_exponent == -2)
  {
    min_float = number;
    max_float = number;
    range_exponent = -1;
    return;
  }
  else
  {
    if (number < min_float)
    {
      min_float = number;
    }
    else if (number > max_float)
    {
      max_float = number;
    }
    else
    {
      return;
    }
  }
  float range = max_float - min_float;
  int exponent = (((U32&)range) >> 23) - ((((U32&)range) & 0x007FFFFF) == 0);
  if (exponent > range_exponent)
  {
    range_exponent = exponent;
    int i,j;
    // the mantissa of numbers with exponent max_exponent will need to
    // preserve the (bits - (range_exponent - max_exponent)) highest
    // bits
//    fprintf(stderr,"exponent range %d\n",range_exponent); 
    i = 22 - bits + range_exponent;
    if (i > 255) i = 255;
    j = 22;
    while (j >= 0 && i >= 0)
    {
      if (rmMantissaNoneLow[i]) {delete rmMantissaNoneLow[i]; rmMantissaNoneLow[i] = 0;};
      if (rmMantissaNoneHigh[i]) {delete rmMantissaNoneHigh[i]; rmMantissaNoneHigh[i] = 0;};
      if (rmMantissaLastLow[i]) {delete rmMantissaLastLow[i]; rmMantissaLastLow[i] = 0;};
      if (rmMantissaLastHigh[i]) {delete rmMantissaLastHigh[i]; rmMantissaLastHigh[i] = 0;};
      if (rmMantissaAcrossLow[i]) {delete rmMantissaAcrossLow[i]; rmMantissaAcrossLow[i] = 0;};
      if (rmMantissaAcrossHigh[i]) {delete rmMantissaAcrossHigh[i]; rmMantissaAcrossHigh[i] = 0;};
      bits_for_mantissa[i] = j;
      i--;
      j--;
    }
    while (i >= 0 && (bits_for_mantissa[i] != -1))
    {
      if (rmMantissaNoneLow[i]) {delete rmMantissaNoneLow[i]; rmMantissaNoneLow[i] = 0;};
      if (rmMantissaNoneHigh[i]) {delete rmMantissaNoneHigh[i]; rmMantissaNoneHigh[i] = 0;};
      if (rmMantissaLastLow[i]) {delete rmMantissaLastLow[i]; rmMantissaLastLow[i] = 0;};
      if (rmMantissaLastHigh[i]) {delete rmMantissaLastHigh[i]; rmMantissaLastHigh[i] = 0;};
      if (rmMantissaAcrossLow[i]) {delete rmMantissaAcrossLow[i]; rmMantissaAcrossLow[i] = 0;};
      if (rmMantissaAcrossHigh[i]) {delete rmMantissaAcrossHigh[i]; rmMantissaAcrossHigh[i] = 0;};
      bits_for_mantissa[i] = -1;
      i--;
    }
  }
}

F32 FloatCompressor::CompressNone(F32 fReal)
{
  num_predictions_none++;

  I32 signReal = (((U32&)fReal) & 0x80000000) == 0x80000000;
  I32 exponentReal = (((U32&)fReal) & 0x7F800000) >> 23;
  I32  mantissaReal = (((U32&)fReal) & 0x007FFFFF);

  if (bits_for_mantissa[exponentReal] >= 0)
  {
    compress_exponent(0,exponentReal,ae_exponent_none,rmExponentNone);
    compress_sign(0,signReal,ae_sign_none,rmSignNone);
    mantissaReal = compress_mantissa(exponentReal,0,mantissaReal,ae_mantissa_none,rmMantissaNoneLow,rmMantissaNoneHigh);

    ((U32&)fReal) = (((U32&)fReal) & 0xFF800000) | mantissaReal;

    if (exponentReal == 255)
    {
      fprintf(stderr,"WARNING: exponentReal is 255 -> infinite number\n");
    }
    else
    {
      if (bits) update_precision(fReal);
    }

    return fReal;
  }
  else // the float is quantized to 0.0f
  {
    compress_exponent(0,0,ae_exponent_none,rmExponentNone);

   return 0.0f;
  }
}

F32 FloatCompressor::DecompressNone()
{
  I32 exponentReal = decompress_exponent(0,ad_exponent_none,rmExponentNone);
  if (bits_for_mantissa[exponentReal] >= 0)
  {
    F32 fReal;
    I32 signReal = decompress_sign(0,ad_sign_none,rmSignNone);
    I32 mantissaReal = decompress_mantissa(exponentReal,0,ad_mantissa_none,rmMantissaNoneLow,rmMantissaNoneHigh);
    if (signReal)
    {
      ((U32&)fReal) = 0x80000000 | (exponentReal << 23) | mantissaReal;
    }
    else
    {
      ((U32&)fReal) = (exponentReal << 23) | mantissaReal;
    }
    if (exponentReal == 255)
    {
      fprintf(stderr,"WARNING: exponentReal is 255 -> infinite number\n");
    }
    else
    {
      if (bits) update_precision(fReal);
    }
    return fReal;
  }
  else
  {
    return 0.0f;
  }
}

F32 FloatCompressor::CompressLast(F32 fPred, F32 fReal)
{
  num_predictions_last++;

  I32 signPred = (((U32&)fPred) & 0x80000000) == 0x80000000;
  I32 exponentPred = (((U32&)fPred) & 0x7F800000) >> 23;
  I32 mantissaPred = (((U32&)fPred) & 0x007FFFFF);

  I32 signReal = (((U32&)fReal) & 0x80000000) == 0x80000000;
  I32 exponentReal = (((U32&)fReal) & 0x7F800000) >> 23;
  I32 mantissaReal = (((U32&)fReal) & 0x007FFFFF);

  if (bits_for_mantissa[exponentReal] >= 0)
  {
    // compress exponent
    compress_exponent(exponentPred,exponentReal,ae_exponent_last,rmExponentLast);
    // compress sign
    if (signPred == signReal)
    {
      compress_sign(exponentReal,0,ae_sign_last,rmSignLast); // correct predicted
    }
    else
    {
      compress_sign(exponentReal,1,ae_sign_last,rmSignLast); // incorrect predicted
    }
    // compress mantissa
    if (signPred == signReal) // sign predicted correctly
    {
      if (exponentPred == exponentReal)
      {
        mantissaReal = compress_mantissa(exponentReal,mantissaPred,mantissaReal,ae_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
      }
      else if (exponentPred > exponentReal)
      {
        mantissaReal = compress_mantissa(exponentReal,0x007FFFFF,mantissaReal,ae_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
      }
      else
      {
        mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
      }
    }
    else if (signPred) // predicted negative, but in reality it is positive
    {
      mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
    }
    else // predicted positive, but in reality it is negative
    {
      mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
    }

    ((U32&)fReal) = (((U32&)fReal) & 0xFF800000) | mantissaReal;

    // update precision
    if (exponentReal == 255)
    {
      fprintf(stderr,"WARNING: exponentReal is 255 -> infinite number\n");
    }
    else
    {
      if (bits) update_precision(fReal);
    }
    return fReal;
  }
  else // the float is quantized to 0.0f
  {
    // compress only exponent
    compress_exponent(exponentPred,0,ae_exponent_last,rmExponentLast);
    return 0.0f;
  }
}

F32 FloatCompressor::DecompressLast(F32 fPred)
{
  I32 exponentPred = (((U32&)fPred) & 0x7F800000) >> 23;
  I32 exponentReal = decompress_exponent(exponentPred,ad_exponent_last,rmExponentLast);

  if (bits_for_mantissa[exponentReal] >= 0)
  {
    F32 fReal;
    I32 signPred = (((U32&)fPred) & 0x80000000) == 0x80000000;
    I32 mantissaPred = (((U32&)fPred) & 0x007FFFFF);
    I32 signReal;
    I32 mantissaReal;
    // decompress sign
    if (decompress_sign(exponentReal,ad_sign_last,rmSignLast))
    {
      signReal = !signPred;
    }
    else
    {
      signReal = signPred;
    }
    // decompress mantissa
    if (signPred == signReal) // sign predicted correctly
    {
      if (exponentPred == exponentReal) // exponent predicted correctly
      {
        mantissaReal = decompress_mantissa(exponentReal,mantissaPred,ad_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
      }
      else if (exponentPred > exponentReal) // overshot exponent prediction
      {
        mantissaReal = decompress_mantissa(exponentReal,0x007FFFFF,ad_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
      }
      else // undershot exponent prediction
      {
        mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
      }
    }
    else if (signPred) // predicted negative, but in reality it is positive
    {
      mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
    }
    else // predicted positive, but in reality it is negative
    {
      mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_last,rmMantissaLastLow,rmMantissaLastHigh);
    }
    // put together decompressed float value
    if (signReal)
    {
      ((U32&)fReal) = 0x80000000 | (exponentReal << 23) | mantissaReal;
    }
    else
    {
      ((U32&)fReal) = (exponentReal << 23) | mantissaReal;
    }
    // update precision
    if (exponentReal == 255)
    {
      fprintf(stderr,"WARNING: exponentReal is 255 -> infinite number\n");
    }
    else
    {
      if (bits) update_precision(fReal);
    }
    return fReal;
  }
  return 0.0f; // we should never get here
}

F32 FloatCompressor::CompressAcross(F32 fPred, F32 fReal)
{
  num_predictions_across++;

  I32 signPred = (((U32&)fPred) & 0x80000000) == 0x80000000;
  I32 exponentPred = (((U32&)fPred) & 0x7F800000) >> 23;
  I32 mantissaPred = (((U32&)fPred) & 0x007FFFFF);

  I32 signReal = (((U32&)fReal) & 0x80000000) == 0x80000000;
  I32 exponentReal = (((U32&)fReal) & 0x7F800000) >> 23;
  I32 mantissaReal = (((U32&)fReal) & 0x007FFFFF);

  if (bits_for_mantissa[exponentReal] >= 0)
  {
    // compress exponent
    compress_exponent(exponentPred,exponentReal,ae_exponent_across,rmExponentAcross);
    // compress sign
    if (signPred == signReal)
    {
      compress_sign(exponentReal,0,ae_sign_across,rmSignAcross); // correct predicted
    }
    else
    {
      compress_sign(exponentReal,1,ae_sign_across,rmSignAcross); // incorrect predicted
    }
    // compress mantissa
    if (signPred == signReal) // sign predicted correctly
    {
      if (exponentPred == exponentReal)
      {
        mantissaReal = compress_mantissa(exponentReal,mantissaPred,mantissaReal,ae_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
      }
      else if (exponentPred > exponentReal)
      {
        mantissaReal = compress_mantissa(exponentReal,0x007FFFFF,mantissaReal,ae_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
      }
      else
      {
        mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
      }
    }
    else if (signPred) // predicted negative, but in reality it is positive
    {
      mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
    }
    else // predicted positive, but in reality it is negative
    {
      mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
    }

    ((U32&)fReal) = (((U32&)fReal) & 0xFF800000) | mantissaReal;

    // update precision
    if (exponentReal == 255)
    {
      fprintf(stderr,"WARNING: exponentReal is 255 -> infinite number\n");
    }
    else
    {
      if (bits) update_precision(fReal);
    }

    return fReal;
  }
  else // the float is quantized to 0.0f
  {
    // compress only exponent
    compress_exponent(exponentPred,0,ae_exponent_across,rmExponentAcross);

    return 0.0f;
  }
}

F32 FloatCompressor::DecompressAcross(F32 fPred)
{
  I32 exponentPred = (((U32&)fPred) & 0x7F800000) >> 23;
  I32 exponentReal = decompress_exponent(exponentPred,ad_exponent_across,rmExponentAcross);

  if (bits_for_mantissa[exponentReal] >= 0)
  {
    F32 fReal;
    I32 signPred = (((U32&)fPred) & 0x80000000) == 0x80000000;
    I32 mantissaPred = (((U32&)fPred) & 0x007FFFFF);
    I32 signReal;
    I32 mantissaReal;
    // decompress sign
    if (decompress_sign(exponentReal,ad_sign_across,rmSignAcross))
    {
      signReal = !signPred;
    }
    else
    {
      signReal = signPred;
    }
    // decompress mantissa
    if (signPred == signReal) // sign predicted correctly
    {
      if (exponentPred == exponentReal) // exponent predicted correctly
      {
        mantissaReal = decompress_mantissa(exponentReal,mantissaPred,ad_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
      }
      else if (exponentPred > exponentReal) // overshot exponent prediction
      {
        mantissaReal = decompress_mantissa(exponentReal,0x007FFFFF,ad_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
      }
      else // undershot exponent prediction
      {
        mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
      }
    }
    else if (signPred) // predicted negative, but in reality it is positive
    {
      mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
    }
    else // predicted positive, but in reality it is negative
    {
      mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_across,rmMantissaAcrossLow,rmMantissaAcrossHigh);
    }
    // put together decompressed float value
    if (signReal)
    {
      ((U32&)fReal) = 0x80000000 | (exponentReal << 23) | mantissaReal;
    }
    else
    {
      ((U32&)fReal) = (exponentReal << 23) | mantissaReal;
    }
    // update precision
    if (exponentReal == 255)
    {
      fprintf(stderr,"WARNING: exponentReal is 255 -> infinite number\n");
    }
    else
    {
      if (bits) update_precision(fReal);
    }
    return fReal;
  }
  return 0.0f; // we should never get here
}

F32 FloatCompressor::CompressWithin(F32 fPred, F32 fReal)
{
  num_predictions_within++;

  I32 signPred = (((U32&)fPred) & 0x80000000) == 0x80000000;
  I32 exponentPred = (((U32&)fPred) & 0x7F800000) >> 23;
  I32 mantissaPred = (((U32&)fPred) & 0x007FFFFF);

  I32 signReal = (((U32&)fReal) & 0x80000000) == 0x80000000;
  I32 exponentReal = (((U32&)fReal) & 0x7F800000) >> 23;
  I32 mantissaReal = (((U32&)fReal) & 0x007FFFFF);

  if (bits_for_mantissa[exponentReal] >= 0)
  {
    // compress exponent
    compress_exponent(exponentPred,exponentReal,ae_exponent_within,rmExponentWithin);
    // compress sign
    if (signPred == signReal)
    {
      compress_sign(exponentReal,0,ae_sign_within,rmSignWithin); // correct predicted
    }
    else
    {
      compress_sign(exponentReal,1,ae_sign_within,rmSignWithin); // incorrect predicted
    }
    // compress mantissa
    if (signPred == signReal) // sign predicted correctly
    {
      if (exponentPred == exponentReal)
      {
        mantissaReal = compress_mantissa(exponentReal,mantissaPred,mantissaReal,ae_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
      }
      else if (exponentPred > exponentReal)
      {
        mantissaReal = compress_mantissa(exponentReal,0x007FFFFF,mantissaReal,ae_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
      }
      else
      {
        mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
      }
    }
    else if (signPred) // predicted negative, but in reality it is positive
    {
      mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
    }
    else // predicted positive, but in reality it is negative
    {
      mantissaReal = compress_mantissa(exponentReal,0x0,mantissaReal,ae_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
    }

    ((U32&)fReal) = (((U32&)fReal) & 0xFF800000) | mantissaReal;

    // update precision
    if (exponentReal == 255)
    {
      fprintf(stderr,"WARNING: exponentReal is 255 -> infinite number\n");
    }
    else
    {
      if (bits) update_precision(fReal);
    }
    return fReal;
  }
  else // the float is quantized to 0.0f
  {
    // compress only exponent
    compress_exponent(exponentPred,0,ae_exponent_within,rmExponentWithin);
    return 0.0f;
  }
}

F32 FloatCompressor::DecompressWithin(F32 fPred)
{
  I32 exponentPred = (((U32&)fPred) & 0x7F800000) >> 23;
  I32 exponentReal = decompress_exponent(exponentPred,ad_exponent_within,rmExponentWithin);

  if (bits_for_mantissa[exponentReal] >= 0)
  {
    F32 fReal;
    I32 signPred = (((U32&)fPred) & 0x80000000) == 0x80000000;
    I32 mantissaPred = (((U32&)fPred) & 0x007FFFFF);
    I32 signReal;
    I32 mantissaReal;
    // decompress sign
    if (decompress_sign(exponentReal,ad_sign_within,rmSignWithin))
    {
      signReal = !signPred;
    }
    else
    {
      signReal = signPred;
    }
    // decompress mantissa
    if (signPred == signReal) // sign predicted correctly
    {
      if (exponentPred == exponentReal) // exponent predicted correctly
      {
        mantissaReal = decompress_mantissa(exponentReal,mantissaPred,ad_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
      }
      else if (exponentPred > exponentReal) // overshot exponent prediction
      {
        mantissaReal = decompress_mantissa(exponentReal,0x007FFFFF,ad_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
      }
      else // undershot exponent prediction
      {
        mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
      }
    }
    else if (signPred) // predicted negative, but in reality it is positive
    {
      mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
    }
    else // predicted positive, but in reality it is negative
    {
      mantissaReal = decompress_mantissa(exponentReal,0x0,ad_mantissa_within,rmMantissaWithinLow,rmMantissaWithinHigh);
    }
    // put together decompressed float value
    if (signReal)
    {
      ((U32&)fReal) = 0x80000000 | (exponentReal << 23) | mantissaReal;
    }
    else
    {
      ((U32&)fReal) = (exponentReal << 23) | mantissaReal;
    }
    // update precision
    if (exponentReal == 255)
    {
      fprintf(stderr,"WARNING: exponentReal is 255 -> infinite number\n");
    }
    else
    {
      if (bits) update_precision(fReal);
    }
    return fReal;
  }
  return 0.0f; // we should never get here
}
