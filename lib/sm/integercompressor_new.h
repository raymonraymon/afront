/*
===============================================================================

  FILE:  integercompressor_new.h
  
  CONTENTS:
 
    This compressor provides three different contexts for encoding integer
    numbers whose range is confined to lie between 2 and 24 bits, which is
    specified with the SetPrecision function. Two of the encoding functions
    take a integer prediction as input. The other will predict the integer
    using the last integer that was encoded.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    09 January 2005 -- completed the bit table of setup_bits()
    27 July 2004 -- the higher order bits should get the bigger tables
    08 January 2004 -- created after clarifying the travel reimbursement claim
  
===============================================================================
*/
#ifndef INTEGER_COMPRESSOR_NEW_H
#define INTEGER_COMPRESSOR_NEW_H

#include "mydefs.h"

#include "rangeencoder.h"
#include "rangedecoder.h"
#include "rangemodel.h"

class IntegerCompressorNew
{
public:

  // SetPrecision:
  void SetPrecision(I32 iBits);
  // GetPrecision:
  I32 GetPrecision();

  // SetRange:
  void SetRange(I32 iRange);
  // GetRange:
  I32 GetRange();

  // SetupCompressor:
  void SetupCompressor(RangeEncoder* re);
  void FinishCompressor();

  // Compress:
  void CompressNone(I32 iReal);
  void CompressLast(I32 iPred, I32 iReal);
  void CompressAcross(I32 iPred, I32 iReal);

  // SetupDecompressor:
  void SetupDecompressor(RangeDecoder* rd);
  void FinishDecompressor();

  // Deompress:
  I32 DecompressNone();
  I32 DecompressLast(I32 iPred);
  I32 DecompressAcross(I32 iPred);

  // Constructor:
  IntegerCompressorNew();
  // Destructor:
  ~IntegerCompressorNew();

  int num_predictions_small; // for statistics only

private:
  // Private Functions

  void setup_bits();
  void writeCorrector(I32 corr, RangeEncoder* ae, RangeModel* amSmall, RangeModel* amHigh);
  I32 readCorrector(RangeDecoder* ad, RangeModel* amSmall, RangeModel* amHigh);

  // Private Variables
  int bits;
  int range;

  int corr_range;
  int corr_max;
  int corr_min;

  int bitsSmall;
  int bitsHigh;
  int bitsLow;

  int smallCutoff;
  int lowMask;
  int highPosOffset;
  int highNegOffset;
  int highestNegative;

  int last;

  RangeEncoder* ae_none;
  RangeEncoder* ae_last;
  RangeEncoder* ae_across;

  RangeDecoder* ad_none;
  RangeDecoder* ad_last;
  RangeDecoder* ad_across;

  RangeModel* amLowPos;
  RangeModel* amLowNeg;

  RangeModel* amSmallNone;
  RangeModel* amHighNone;

  RangeModel* amSmallLast;
  RangeModel* amHighLast;

  RangeModel* amSmallAcross;
  RangeModel* amHighAcross;
};

#endif
