/*
===============================================================================

  FILE:  integercompressor_new.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see corresponding header file
  
===============================================================================
*/

#include "integercompressor_new.h"

IntegerCompressorNew::IntegerCompressorNew()
{
  bits = 16;
  range = 0;

  bitsSmall = 9;
  bitsHigh = 12;
  bitsLow = 4;

  smallCutoff = 0;
  lowMask = 0;
  highPosOffset = 0;
  highNegOffset = 0;
  highestNegative = 0;

  last = 0;

  num_predictions_small = 0;

  ae_none = 0;
  ae_last = 0;
  ae_across = 0;

  ad_none = 0;
  ad_last = 0;
  ad_across = 0;

  amLowPos = 0;
  amLowNeg = 0;

  amSmallNone = 0;
  amHighNone = 0;

  amSmallLast = 0;
  amHighLast = 0;

  amSmallAcross = 0;
  amHighAcross = 0;
}

IntegerCompressorNew::~IntegerCompressorNew()
{
}

void IntegerCompressorNew::setup_bits()
{
  if (bits == 24)
  {
    bitsSmall = 14;
    bitsHigh = 12;
    bitsLow = 12;
  }
  else if (bits == 23)
  {
    bitsSmall = 14;
    bitsHigh = 12;
    bitsLow = 11;
  }
  else if (bits == 22)
  {
    bitsSmall = 13;
    bitsHigh = 12;
    bitsLow = 10;
  }
  else if (bits == 21)
  {
    bitsSmall = 12;
    bitsHigh = 12;
    bitsLow = 9;
  }
  else if (bits == 20)
  {
    bitsSmall = 11;
    bitsHigh = 12;
    bitsLow = 8;
  }
  else if (bits == 19)
  {
    bitsSmall = 10;
    bitsHigh = 12;
    bitsLow = 7;
  }
  else if (bits == 18)
  {
    bitsSmall = 10;
    bitsHigh = 12;
    bitsLow = 6;
  }
  else if (bits == 17)
  {
    bitsSmall = 9;
    bitsHigh = 12;
    bitsLow = 5;
  }
  else if (bits == 16)
  {
    bitsSmall = 9;
    bitsHigh = 12;
    bitsLow = 4;
  }
  else if (bits == 15)
  {
    bitsSmall = 8;
    bitsHigh = 12;
    bitsLow = 3;
  }
  else if (bits == 14)
  {
    bitsSmall = 8;
    bitsHigh = 14;
    bitsLow = 0;
  }
  else if (bits == 13)
  {
    bitsSmall = 7;
    bitsHigh = 13;
    bitsLow = 0;
  }
  else if (bits == 12)
  {
    bitsSmall = 7;
    bitsHigh = 12;
    bitsLow = 0;
  }
  else if (bits == 11)
  {
    bitsSmall = 6;
    bitsHigh = 11;
    bitsLow = 0;
  }
  else if (bits == 10)
  {
    bitsSmall = 5;
    bitsHigh = 10;
    bitsLow = 0;
  }
  else if (bits == 9)
  {
    bitsSmall = 5;
    bitsHigh = 9;
    bitsLow = 0;
  }
  else if (bits == 8)
  {
    bitsSmall = 5;
    bitsHigh = 8;
    bitsLow = 0;
  }
  else if (bits == 7)
  {
    bitsSmall = 4;
    bitsHigh = 7;
    bitsLow = 0;
  }
  else if (bits < 7)
  {
    bitsSmall = bits;
    bitsHigh = 0;
    bitsLow = 0;
  }
}

void IntegerCompressorNew::SetupCompressor(RangeEncoder* re)
{
  ae_none = re;
  ae_last = re;
  ae_across = re;

  setup_bits();

  // the correctors must fall into this interval

  if (range)
  {
    corr_range = range;
  }
  else
  {
    corr_range = 1 << bits;
  }
  corr_min = -(corr_range/2);
  corr_max = corr_min + corr_range - 1;

  if (corr_range <= (1 << bitsSmall))
  {
    // sometimes (especially for low precision) the entire corr_range falls under the small cutoff

    if (bitsSmall > 0)
    {
      smallCutoff = 1 + corr_range/2;
    }
    else
    {
      smallCutoff = 0;
    }

    amSmallNone = new RangeModel(1 + corr_range,0,1,20000,16); 
    amSmallLast = new RangeModel(1 + corr_range,0,1,20000,16); 
    amSmallAcross = new RangeModel(1 + corr_range,0,1,20000,16); 
  }
  else
  {
    // usually we hope that still a large number of correctors fall under the small cutoff

    if (bitsSmall > 0)
    {
      smallCutoff = (1 << (bitsSmall-1));
    }
    else
    {
      smallCutoff = 0;
    }

    amSmallNone = new RangeModel((1 << bitsSmall),0,1,20000,16); 
    amSmallLast = new RangeModel((1 << bitsSmall),0,1,20000,16); 
    amSmallAcross = new RangeModel((1 << bitsSmall),0,1,20000,16); 

    // only if they do not we need the other stuff

    if (bitsLow == 0)
    {
      highNegOffset = -(corr_min);
      highPosOffset = -(corr_min + (2*smallCutoff - 1));

      amHighNone = new RangeModel(corr_range - (2*smallCutoff - 1),0,1,20000,16);
      amHighLast = new RangeModel(corr_range - (2*smallCutoff - 1),0,1,20000,16);
      amHighAcross = new RangeModel(corr_range - (2*smallCutoff - 1),0,1,20000,16);
    }
    else
    {
      lowMask = (1 << bitsLow) - 1;

      amLowPos = new RangeModel((1 << bitsLow),0,1,20000,16);
      amLowNeg = new RangeModel((1 << bitsLow),0,1,20000,16);

      highNegOffset = -(corr_min >> bitsLow);
      highPosOffset = -(corr_min >> bitsLow);

      amHighNone = new RangeModel((corr_max >> bitsLow) - (corr_min >> bitsLow) + 1,0,1,20000,16);
      amHighLast = new RangeModel((corr_max >> bitsLow) - (corr_min >> bitsLow) + 1,0,1,20000,16);
      amHighAcross = new RangeModel((corr_max >> bitsLow) - (corr_min >> bitsLow) + 1,0,1,20000,16);
    }
  }
}

void IntegerCompressorNew::FinishCompressor()
{
  if (amLowPos) delete amLowPos;
  if (amLowNeg) delete amLowNeg;

  if (amSmallNone) delete amSmallNone;
  if (amHighNone) delete amHighNone;

  if (amSmallLast) delete amSmallLast;
  if (amHighLast) delete amHighLast;

  if (amSmallAcross) delete amSmallAcross;
  if (amHighAcross) delete amHighAcross;
}

void IntegerCompressorNew::SetupDecompressor(RangeDecoder* rd)
{
  ad_none = rd;
  ad_last = rd;
  ad_across = rd;

  setup_bits();

  // the correctors must fall into this interval

  if (range)
  {
    corr_range = range;
  }
  else
  {
    corr_range = 1 << bits;
  }
  corr_min = -(corr_range/2);
  corr_max = corr_min + corr_range - 1;

  if (corr_range <= (1 << bitsSmall))
  {
    // sometimes (especially for low precision) the entire corr_range falls under the small cutoff

    smallCutoff = 1 + corr_range/2;

    amSmallNone = new RangeModel(1 + corr_range,0,0,20000,16); 
    amSmallLast = new RangeModel(1 + corr_range,0,0,20000,16); 
    amSmallAcross = new RangeModel(1 + corr_range,0,0,20000,16); 
  }
  else
  {
    // usually we hope that still a large number of correctors fall under the small cutoff

    smallCutoff = (1 << (bitsSmall-1));

    amSmallNone = new RangeModel((1 << bitsSmall),0,0,20000,16); 
    amSmallLast = new RangeModel((1 << bitsSmall),0,0,20000,16); 
    amSmallAcross = new RangeModel((1 << bitsSmall),0,0,20000,16); 

    // only if they do not we need the other stuff

    if (bitsLow == 0)
    {
      highNegOffset = -(corr_min);
      highPosOffset = -(corr_min + (2*smallCutoff - 1));
      highestNegative = -(corr_min + smallCutoff);

      amHighNone = new RangeModel(corr_range - (2*smallCutoff - 1),0,0,20000,16);
      amHighLast = new RangeModel(corr_range - (2*smallCutoff - 1),0,0,20000,16);
      amHighAcross = new RangeModel(corr_range - (2*smallCutoff - 1),0,0,20000,16);
    }
    else
    {
      lowMask = (1 << bitsLow) - 1;

      amLowPos = new RangeModel((1 << bitsLow),0,0,20000,16);
      amLowNeg = new RangeModel((1 << bitsLow),0,0,20000,16);

      highNegOffset = -(corr_min >> bitsLow);
      highPosOffset = -(corr_min >> bitsLow);
      highestNegative = -(corr_min >> bitsLow);

      amHighNone = new RangeModel((corr_max >> bitsLow) - (corr_min >> bitsLow) + 1,0,0,20000,16);
      amHighLast = new RangeModel((corr_max >> bitsLow) - (corr_min >> bitsLow) + 1,0,0,20000,16);
      amHighAcross = new RangeModel((corr_max >> bitsLow) - (corr_min >> bitsLow) + 1,0,0,20000,16);
    }
  }
}

void IntegerCompressorNew::FinishDecompressor()
{
  if (amLowPos) delete amLowPos;
  if (amLowNeg) delete amLowNeg;

  if (amSmallNone) delete amSmallNone;
  if (amHighNone) delete amHighNone;

  if (amSmallLast) delete amSmallLast;
  if (amHighLast) delete amHighLast;

  if (amSmallAcross) delete amSmallAcross;
  if (amHighAcross) delete amHighAcross;
}

//-----------------------------------------------------------------------------
// SetPrecision:
//-----------------------------------------------------------------------------
void IntegerCompressorNew::SetPrecision(I32 iBits)
{
  bits = iBits;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SetRange:
//-----------------------------------------------------------------------------
void IntegerCompressorNew::SetRange(I32 iRange)
{
  range = iRange;
}

void IntegerCompressorNew::writeCorrector(I32 corr, RangeEncoder* ae, RangeModel* amSmall, RangeModel* amHigh)
{
  if (corr < 0)
  {
    if (-corr < smallCutoff)
    {
      num_predictions_small++;
      ae->encode(amSmall,corr + smallCutoff);
    }
    else
    {
      ae->encode(amSmall, 0); // does not fit into the small
      if (bitsLow == 0)
      {
        ae->encode(amHigh, corr + highNegOffset);
      }
      else
      {
        ae->encode(amHigh, (corr >> bitsLow) + highNegOffset);
        ae->encode(amLowNeg, corr & lowMask);
      }
    }
  }
  else
  {
    if (corr < smallCutoff)
    {
      num_predictions_small++;
      ae->encode(amSmall,corr + smallCutoff);
    }
    else
    {
      ae->encode(amSmall, 0); // does not fit into the small
      if (bitsLow == 0)
      {
        ae->encode(amHigh, corr + highPosOffset);
      }
      else
      {
        ae->encode(amHigh, (corr >> bitsLow) + highPosOffset);
        ae->encode(amLowPos, corr & lowMask);
      }
    }
  }
}

void IntegerCompressorNew::CompressNone(I32 real)
{
  int corr = real - last;
  if (corr < corr_min) corr += corr_range;
  else if (corr > corr_max) corr -= corr_range;
  writeCorrector(corr, ae_none, amSmallNone, amHighNone);
  last = real;
}

void IntegerCompressorNew::CompressLast(I32 pred, I32 real)
{
  int corr = real - pred;
  if (corr < corr_min) corr += corr_range;
  else if (corr > corr_max) corr -= corr_range;
  writeCorrector(corr, ae_last, amSmallLast, amHighLast);
  last = real;
}

void IntegerCompressorNew::CompressAcross(I32 pred, I32 real)
{
  int corr = real - pred;
  if (corr < corr_min) corr += corr_range;
  else if (corr > corr_max) corr -= corr_range;
  writeCorrector(corr, ae_across, amSmallAcross, amHighAcross);
  last = real;
}

I32 IntegerCompressorNew::readCorrector(RangeDecoder* ad, RangeModel* amSmall, RangeModel* amHigh)
{
  int corr = ad->decode(amSmall);

  if (corr)
  {
    return corr - smallCutoff;
  }
  else
  {
    if (bitsLow == 0)
    {
      corr = ad->decode(amHigh);
      if (corr <= highestNegative)
      {
        return corr - highNegOffset;
      }
      else
      {
        return corr - highPosOffset;
      }
    }
    else
    {
      corr = ad->decode(amHigh);
      if (corr <= highestNegative)
      {
        corr = (corr - highNegOffset) << bitsLow;
        return corr | ad->decode(amLowNeg);
      }
      else
      {
        corr = (corr - highPosOffset) << bitsLow;
        return corr | ad->decode(amLowPos);
      }
    }
  }
}

I32 IntegerCompressorNew::DecompressNone()
{
  int real = last + readCorrector(ad_none, amSmallNone, amHighNone);
  if (real < 0) real += corr_range;
  else if (real >= corr_range) real -= corr_range;
  last = real;
  return real;
}

I32 IntegerCompressorNew::DecompressLast(I32 pred)
{
  int real = pred + readCorrector(ad_last, amSmallLast, amHighLast);
  if (real < 0) real += corr_range;
  else if (real >= corr_range) real -= corr_range;
  last = real;
  return real;
}

I32 IntegerCompressorNew::DecompressAcross(I32 pred)
{
  int real = pred + readCorrector(ad_across, amSmallAcross, amHighAcross);
  if (real < 0) real += corr_range;
  else if (real >= corr_range) real -= corr_range;
  last = real;
  return real;
}
