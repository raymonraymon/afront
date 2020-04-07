/*
===============================================================================

  FILE:  integercompressor.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    Copyright (C) 2000  Martin Isenburg (isenburg@cs.unc.edu)
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see corresponding header file
  
===============================================================================
*/

#include "integercompressor.h"

IntegerCompressor::IntegerCompressor()
{
  bits = 16;

  bitsBigLow = 0;
  bitsBigHigh = 0;
  bitsSmall = 0;

  bitsMaskBigLow = 0;
  bitsMaskSmall = 0;

  maxNone = 0;
  maxLast = 0;
  maxAcross = 0;

  ae_none = 0;
  ae_last = 0;
  ae_across = 0;

  ad_none = 0;
  ad_last = 0;
  ad_across = 0;

  amBitsLastBigLow = 0;
  amBitsLastBigHigh = 0;
  amBitsLastSmall = 0;

  amBitsAcrossBigLow = 0;
  amBitsAcrossBigLow = 0;
  amBitsAcrossSmall = 0;

  num_predictions_none = 0;
  num_predictions_last = 0;
  num_predictions_across = 0;
}

IntegerCompressor::~IntegerCompressor()
{
  if (ae_none != ae_last)
  {
    if (ae_none) delete ae_none;
    if (ae_last) delete ae_last;
    if (ae_across) delete ae_across;
  }

  if (ad_none != ad_last)
  {
    if (ad_none) delete ad_none;
    if (ad_last) delete ad_last;
    if (ad_across) delete ad_across;
  }

  if (amBitsLastBigLow) delete amBitsLastBigLow;
  if (amBitsLastBigHigh) delete amBitsLastBigHigh;
  if (amBitsLastSmall) delete amBitsLastSmall;

  if (amBitsAcrossBigLow) delete amBitsAcrossBigLow;
  if (amBitsAcrossBigHigh) delete amBitsAcrossBigHigh;
  if (amBitsAcrossSmall) delete amBitsAcrossSmall;
}

void IntegerCompressor::SetupCompressor(RangeEncoder* re)
{
  if (bits == 24)
  {
    // 14-9 = 23+1 bits
    bitsBigLow = 14;
    bitsBigHigh = 9;
    bitsSmall = 12;
  }
  else if (bits == 22)
  {
    // 13-8 = 21+1 bits
    bitsBigLow = 13;
    bitsBigHigh = 8;
    bitsSmall = 11;
  }
  else if (bits == 20)
  {
    // 13-6 = 19+1 bits
    bitsBigLow = 13;
    bitsBigHigh = 6;
    bitsSmall = 10;
  }
  else if (bits == 18)
  {
    // 12-5 = 17+1 bits
    bitsBigLow = 12;
    bitsBigHigh = 5;
    bitsSmall = 9;
  }
  else if (bits == 16)
  {
    // 11-4 = 15+1 bits
    bitsBigLow = 11;
    bitsBigHigh = 4;
    bitsSmall = 8;
  }
  else if (bits == 14)
  {
    // 10-3 = 13+1 bits
    bitsBigLow = 10;
    bitsBigHigh = 3;
    bitsSmall = 7;
  }
  else if (bits == 12)
  {
    // 11 = 11+1 bits
    bitsBigLow = 11;
    bitsBigHigh = 0;
    bitsSmall = 6;
  }
  else if (bits == 10)
  {
    // 9 = 9+1 bits
    bitsBigLow = 9;
    bitsBigHigh = 0;
    bitsSmall = 5;
  }
  else if (bits == 9)
  {
    // 8 = 8+1 bits
    bitsBigLow = 8;
    bitsBigHigh = 0;
    bitsSmall = 5;
  }
  else if (bits == 8)
  {
    // 4 = 7 bits
    bitsBigLow = 7;
    bitsBigHigh = 0;
    bitsSmall = 4;
  }

  bitsMaskBigLow = (1 << bitsBigLow) - 1;
  bitsMaskSmall = (1 << bitsSmall) - 1;

  amBitsLastBigLow = new RangeModel((1 << bitsBigLow),0,1);
  amBitsAcrossBigLow = new RangeModel((1 << bitsBigLow),0,1);

  if (bitsBigHigh)
  {
    if (((maxLast >> bitsBigLow) + 1) >= (1<<bitsBigHigh))
    {
      amBitsLastBigHigh = new RangeModel((1<<bitsBigHigh),0,1);
    }
    else
    {
      amBitsLastBigHigh = new RangeModel(((maxLast >> bitsBigLow) + 1),0,1);
    }
    if (((maxAcross >> bitsBigLow) + 1) >= (1<<bitsBigHigh))
    {
      amBitsAcrossBigHigh = new RangeModel((1<<bitsBigHigh),0,1);
    }
    else
    {
      amBitsAcrossBigHigh = new RangeModel(((maxAcross >> bitsBigLow) + 1),0,1);
    }
  }
  else
  {
    amBitsLastBigHigh = 0;
    amBitsAcrossBigHigh = 0;
  }

  if ((maxLast + 1) > (1 << bitsSmall))
  {
    amBitsLastSmall = new RangeModel((1 << bitsSmall),0,1);
  }
  else
  {
    amBitsLastSmall = new RangeModel((maxLast + 1),0,1);
  }

  if ((maxAcross + 1) > (1 << bitsSmall))
  {
    amBitsAcrossSmall = new RangeModel((1 << bitsSmall),0,1);
  }
  else
  {
    amBitsAcrossSmall = new RangeModel((maxAcross + 1),0,1);
  }

  if (re == 0)
  {
    ae_none = new RangeEncoder(0);
    ae_last = new RangeEncoder(0);
    ae_across = new RangeEncoder(0);
  }
  else
  {
    ae_none = re;
    ae_last = re;
    ae_across = re;
  }
}

void IntegerCompressor::FinishCompressor()
{
  if (ae_none != ae_last)
  {
    fprintf(stderr,"none %d last %d across %d\n", num_predictions_none, num_predictions_last, num_predictions_across);
    ae_none->done();
    ae_last->done();
    ae_across->done();
    fprintf(stderr,"none:   %d bytes %5.2f bpf\n",ae_none->getNumberChars(),(float)(ae_none->getNumberBits())/(float)(num_predictions_none));
    fprintf(stderr,"last:   %d bytes %5.2f bpf\n",ae_last->getNumberChars(),(float)(ae_last->getNumberBits())/(float)(num_predictions_last));
    fprintf(stderr,"across: %d bytes %5.2f bpf\n",ae_across->getNumberChars(),(float)(ae_across->getNumberBits())/(float)(num_predictions_across));
    fprintf(stderr,"TOTAL:  %d bytes %5.2f bpf\n",ae_none->getNumberChars()+ae_last->getNumberChars()+ae_across->getNumberChars(),(float)(ae_none->getNumberBits()+ae_last->getNumberBits()+ae_across->getNumberBits())/(float)(num_predictions_none+num_predictions_last+num_predictions_across));
  }
}

void IntegerCompressor::SetupDecompressor(RangeDecoder* rd)
{
  if (bits == 24)
  {
    // 14-9 = 23+1 bits
    bitsBigLow = 14;
    bitsBigHigh = 9;
    bitsSmall = 12;
  }
  else if (bits == 22)
  {
    // 13-8 = 21+1 bits
    bitsBigLow = 13;
    bitsBigHigh = 8;
    bitsSmall = 11;
  }
  else if (bits == 20)
  {
    // 13-6 = 19+1 bits
    bitsBigLow = 13;
    bitsBigHigh = 6;
    bitsSmall = 10;
  }
  else if (bits == 18)
  {
    // 12-5 = 17+1 bits
    bitsBigLow = 12;
    bitsBigHigh = 5;
    bitsSmall = 9;
  }
  else if (bits == 16)
  {
    // 11-4 = 15+1 bits
    bitsBigLow = 11;
    bitsBigHigh = 4;
    bitsSmall = 8;
  }
  else if (bits == 14)
  {
    // 10-3 = 13+1 bits
    bitsBigLow = 10;
    bitsBigHigh = 3;
    bitsSmall = 7;
  }
  else if (bits == 12)
  {
    // 11 = 11+1 bits
    bitsBigLow = 11;
    bitsBigHigh = 0;
    bitsSmall = 6;
  }
  else if (bits == 10)
  {
    // 9 = 9+1 bits
    bitsBigLow = 9;
    bitsBigHigh = 0;
    bitsSmall = 5;
  }
  else if (bits == 9)
  {
    // 8 = 8+1 bits
    bitsBigLow = 8;
    bitsBigHigh = 0;
    bitsSmall = 5;
  }
  else if (bits == 8)
  {
    // 4 = 7 bits
    bitsBigLow = 7;
    bitsBigHigh = 0;
    bitsSmall = 4;
  }

  bitsMaskBigLow = (1 << bitsBigLow) - 1;
  bitsMaskSmall = (1 << bitsSmall) - 1;

  amBitsLastBigLow = new RangeModel((1 << bitsBigLow),0,0);
  amBitsAcrossBigLow = new RangeModel((1 << bitsBigLow),0,0);

  if (bitsBigHigh)
  {
    if (((maxLast >> bitsBigLow) + 1) >= (1<<bitsBigHigh))
    {
      amBitsLastBigHigh = new RangeModel((1<<bitsBigHigh),0,0);
    }
    else
    {
      amBitsLastBigHigh = new RangeModel(((maxLast >> bitsBigLow) + 1),0,0);
    }
    if (((maxAcross >> bitsBigLow) + 1) >= (1<<bitsBigHigh))
    {
      amBitsAcrossBigHigh = new RangeModel((1<<bitsBigHigh),0,0);
    }
    else
    {
      amBitsAcrossBigHigh = new RangeModel(((maxAcross >> bitsBigLow) + 1),0,0);
    }
  }
  else
  {
    amBitsLastBigHigh = 0;
    amBitsAcrossBigHigh = 0;
  }

  if ((maxLast + 1) > (1 << bitsSmall))
  {
    amBitsLastSmall = new RangeModel((1 << bitsSmall),0,0);
  }
  else
  {
    amBitsLastSmall = new RangeModel((maxLast + 1),0,0);
  }

  if ((maxAcross + 1) > (1 << bitsSmall))
  {
    amBitsAcrossSmall = new RangeModel((1 << bitsSmall),0,0);
  }
  else
  {
    amBitsAcrossSmall = new RangeModel((maxAcross + 1),0,0);
  }

  if (rd == 0)
  {
    ad_none = new RangeDecoder(ae_none->getChars(),ae_none->getNumberChars());
    ad_last = new RangeDecoder(ae_last->getChars(),ae_last->getNumberChars());
    ad_across = new RangeDecoder(ae_across->getChars(),ae_across->getNumberChars());
  }
  else
  {
    ad_none = rd;
    ad_last = rd;
    ad_across = rd;
  }
}

void IntegerCompressor::FinishDecompressor()
{
  if (ad_none != ad_last)
  {
    ad_none->done();
    ad_last->done();
    ad_across->done();
  }
}

//-----------------------------------------------------------------------------
// SetPrecision:
//-----------------------------------------------------------------------------
void IntegerCompressor::SetPrecision(I32 iBits)
{
  bits = iBits;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// UpdateMax:
//-----------------------------------------------------------------------------
void IntegerCompressor::UpdateMax(I32 iNumNone, I32 iNumLast, I32 iNumAcross)
{
  if (iNumNone > maxNone) maxNone = iNumNone;
  if (iNumLast > maxLast) maxLast = iNumLast;
  if (iNumAcross > maxAcross) maxAcross = iNumAcross;
}

//-----------------------------------------------------------------------------
// SetMinMax:
//-----------------------------------------------------------------------------
void IntegerCompressor::SetMax(I32 iMaxNone, I32 iMaxLast, I32 iMaxAcross)
{
  maxNone = iMaxNone;
  maxLast = iMaxLast;
  maxAcross = iMaxAcross;
}

void IntegerCompressor::writeCorrector(I32 corr, RangeEncoder* ae, RangeModel* amBitsSmall, RangeModel* amBitsBigLow, RangeModel* amBitsBigHigh)
{
  int c;
  int sign, low;

  c = corr;
  if (c > 0)      // positive
  {
    sign = 0;
  }
  else if (c < 0)  // negative
  {
    sign = 1;
    c = -c;
  }
  else            // zero
  {
    sign = -1;
  }
  if (c >= bitsMaskSmall)
  {
    ae->encode(amBitsSmall,bitsMaskSmall);
    if (amBitsBigHigh)
    {
      low = c & bitsMaskBigLow;
      c >>= bitsBigLow;
      ae->encode(amBitsBigHigh, c);
      ae->encode(amBitsBigLow, low);
    }
    else
    {
      ae->encode(amBitsBigLow, c);
    }
  }
  else
  {
    ae->encode(amBitsSmall,c);
  }
  if (sign != -1) // positive or negative (e.g. not null)
  {
    ae->encode(2, sign);
  }
}

void IntegerCompressor::CompressNone(I32 pos)
{
  num_predictions_none++;
  ae_none->encode(maxNone, pos);
}

void IntegerCompressor::CompressLast(I32 corr)
{
  num_predictions_last++;
  writeCorrector(corr, ae_last, amBitsLastSmall, amBitsLastBigLow, amBitsLastBigHigh);
}

void IntegerCompressor::CompressAcross(I32 corr)
{
  num_predictions_across++;
  writeCorrector(corr, ae_across, amBitsAcrossSmall, amBitsAcrossBigLow, amBitsAcrossBigHigh);
}

I32 IntegerCompressor::readCorrector(RangeDecoder* ad, RangeModel* amBitsSmall, RangeModel* amBitsBigLow, RangeModel* amBitsBigHigh)
{
  int c;

  c = ad->decode(amBitsSmall);
  if (c == bitsMaskSmall)
  {
    if (amBitsBigHigh)
    {
      c = ad->decode(amBitsBigHigh) << bitsBigLow;
      c |= ad->decode(amBitsBigLow);
    }
    else
    {
      c = ad->decode(amBitsBigLow);
    }
  }
  if (c != 0) // decode sign
  {
    if (ad->decode(2) == 1)
    {
      c = -c;
    }
  }
  return c;
}

I32 IntegerCompressor::DecompressNone()
{
  return ad_none->decode(maxNone);
}

I32 IntegerCompressor::DecompressLast()
{
  return readCorrector(ad_last, amBitsLastSmall, amBitsLastBigLow, amBitsLastBigHigh);
}

I32 IntegerCompressor::DecompressAcross()
{
  return readCorrector(ad_across, amBitsAcrossSmall, amBitsAcrossBigLow, amBitsAcrossBigHigh);
}

