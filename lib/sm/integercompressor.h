/*
===============================================================================

  FILE:  integercompressor.h
  
  CONTENTS:
  
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2004  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    8 January 2004 -- created after clarifying the travel reimbursement claim
  
===============================================================================
*/
#ifndef INTEGER_COMPRESSOR_H
#define INTEGER_COMPRESSOR_H

#include "mydefs.h"

#include "rangeencoder.h"
#include "rangedecoder.h"
#include "rangemodel.h"

class IntegerCompressor
{
public:

  // SetPrecision:
  void SetPrecision(I32 iBits);
  // GetPrecision:
  I32 GetPrecision();

  // UpdateMax:
  void UpdateMax(I32 iNumNone, I32 iNumLast, I32 iNumAcross);
  // SetMax:
  void SetMax(I32 iNumNone, I32 iNumLast, I32 iNumAcross);

  // SetupCompressor:
  void SetupCompressor(RangeEncoder* re);
  void FinishCompressor();

  // Compress:
  void CompressNone(I32 iNumNone);
  void CompressLast(I32 iNumLast);
  void CompressAcross(I32 iNumAcross);

  // SetupDecompressor:
  void SetupDecompressor(RangeDecoder* rd);
  void FinishDecompressor();

  // Deompress:
  I32 DecompressNone();
  I32 DecompressLast();
  I32 DecompressAcross();

  // Constructor:
  IntegerCompressor();
  // Destructor:
  ~IntegerCompressor();

private:
  // Private Functions

  void writeCorrector(I32 corr, RangeEncoder* ae, RangeModel* amBitsSmall, RangeModel* amBitsBigLow, RangeModel* amBitsBigHigh);
  I32 readCorrector(RangeDecoder* ad, RangeModel* amBitsSmall, RangeModel* amBitsBigLow, RangeModel* amBitsBigHigh);

  // Private Variables
  RangeEncoder* ae_none;
  RangeEncoder* ae_last;
  RangeEncoder* ae_across;

  RangeDecoder* ad_none;
  RangeDecoder* ad_last;
  RangeDecoder* ad_across;

  RangeModel* amBitsLastBigLow;
  RangeModel* amBitsLastBigHigh;
  RangeModel* amBitsLastSmall;

  RangeModel* amBitsAcrossBigLow;
  RangeModel* amBitsAcrossBigHigh;
  RangeModel* amBitsAcrossSmall;

  int bits;

  int bitsBigLow;
  int bitsBigHigh;
  int bitsSmall;

  int bitsMaskBigLow;
  int  bitsMaskSmall;

  int maxNone;
  int maxLast;
  int maxAcross;

  int num_predictions_none;
  int num_predictions_last;
  int num_predictions_across;
};

#endif
