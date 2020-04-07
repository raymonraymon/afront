/*
===============================================================================

  FILE:  floatcompressor.h
  
  CONTENTS:
  
    This floating point compressor can:

    * compress floating point numbers in a lossless manner while exploiting
      the benefits of predictive coding (set precision to be 0)

    * uniformly quantize and compress floating point numbers without a-priori
      informtion about the bounding box while guaranteeing k bits of precision
  
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
#ifndef FLOAT_COMPRESSOR_H
#define FLOAT_COMPRESSOR_H

#include "mydefs.h"

#include "rangeencoder.h"
#include "rangedecoder.h"
#include "rangemodel.h"

class FloatCompressor
{
public:

  // SetPrecision:
  inline void SetPrecision(I32 iBits);
  // GetPrecision:
  inline I32 GetPrecision();

  // UpdateMinMax:
  inline void UpdateMinMax(F32 fFloat);
  // SetMinMax:
  inline void SetMinMax(F32 fMin, F32 fMax);
  // GetMinMax:
  inline void GetMinMax(F32& fMin, F32& fMax);

  // SetupCompressor:
  void SetupCompressor(RangeEncoder* re, bool within);
  void FinishCompressor(bool within);

  // Compress:
  F32 CompressNone(F32 fReal);
  F32 CompressLast(F32 fPred, F32 fReal);
  F32 CompressAcross(F32 fPred, F32 fReal);
  F32 CompressWithin(F32 fPred, F32 fReal);

  // Decompress:
  void SetupDecompressor(RangeDecoder* rd, bool within);
  void FinishDecompressor(bool within);

  // Compress:
  F32 DecompressNone();
  F32 DecompressLast(F32 fPred);
  F32 DecompressAcross(F32 fPred);
  F32 DecompressWithin(F32 fPred);

  // Constructor:
  FloatCompressor();
  // Destructor:
  ~FloatCompressor();

private:
  // Private Functions
  void alloc_range_tables(bool within);
  void dealloc_range_tables(bool within);

  void reset_precision();
  void update_precision(float number);

  void compress_exponent(int exponentPred, int exponentReal, RangeEncoder* re_exponent, RangeModel** rmExponent);
  int decompress_exponent(int exponentPred, RangeDecoder* rd_exponent, RangeModel** rmExponent);

  void compress_sign(int exponent, int sign, RangeEncoder* re_sign, RangeModel** rmSign);
  int decompress_sign(int exponent, RangeDecoder* rd_sign, RangeModel** rmSign);

  int compress_mantissa(int exponent, int mantissaPred, int mantissaReal, RangeEncoder* re_mantissa, RangeModel** rmMantissaHigh, RangeModel** rmMantissaLow);
  int decompress_mantissa(int exponent, int mantissaPred, RangeDecoder* rd_mantissa, RangeModel** rmMantissaHigh, RangeModel** rmMantissaLow);

  // Private Variables
  RangeEncoder* ae_sign_none;
  RangeEncoder* ae_exponent_none;
  RangeEncoder* ae_mantissa_none;

  RangeEncoder* ae_sign_last;
  RangeEncoder* ae_exponent_last;
  RangeEncoder* ae_mantissa_last;

  RangeEncoder* ae_sign_across;
  RangeEncoder* ae_exponent_across;
  RangeEncoder* ae_mantissa_across;

  RangeEncoder* ae_sign_within;
  RangeEncoder* ae_exponent_within;
  RangeEncoder* ae_mantissa_within;

  RangeDecoder* ad_sign_none;
  RangeDecoder* ad_exponent_none;
  RangeDecoder* ad_mantissa_none;

  RangeDecoder* ad_sign_last;
  RangeDecoder* ad_exponent_last;
  RangeDecoder* ad_mantissa_last;

  RangeDecoder* ad_sign_across;
  RangeDecoder* ad_exponent_across;
  RangeDecoder* ad_mantissa_across;

  RangeDecoder* ad_sign_within;
  RangeDecoder* ad_exponent_within;
  RangeDecoder* ad_mantissa_within;

  RangeModel** rmSignNone;
  RangeModel** rmExponentNone;
  RangeModel** rmMantissaNoneLow;
  RangeModel** rmMantissaNoneHigh;

  RangeModel** rmSignLast;
  RangeModel** rmExponentLast;
  RangeModel** rmMantissaLastLow;
  RangeModel** rmMantissaLastHigh;

  RangeModel** rmSignAcross;
  RangeModel** rmExponentAcross;
  RangeModel** rmMantissaAcrossLow;
  RangeModel** rmMantissaAcrossHigh;

  RangeModel** rmSignWithin;
  RangeModel** rmExponentWithin;
  RangeModel** rmMantissaWithinLow;
  RangeModel** rmMantissaWithinHigh;

  RangeModel** rmZeroOrExtremeValue;
  RangeModel** rmPositiveOrNegativeValue;

  int* bits_for_mantissa;
  int range_exponent;
  float min_float;
  float max_float;
  int bits;
};

//-----------------------------------------------------------------------------
// SetPrecision:
//-----------------------------------------------------------------------------
inline void FloatCompressor::SetPrecision(I32 iBits)
{
  bits = iBits;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GetPrecision:
//-----------------------------------------------------------------------------
inline I32 FloatCompressor::GetPrecision()
{
  return bits;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// UpdateMinMax:
//-----------------------------------------------------------------------------
inline void FloatCompressor::UpdateMinMax(F32 pPos)
{
  fprintf(stderr,"not yet implemented\n");
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SetMinMax:
//-----------------------------------------------------------------------------
inline void FloatCompressor::SetMinMax(F32 pMin, F32 pMax)
{
  fprintf(stderr,"not yet implemented\n");
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GetMinMax:
//-----------------------------------------------------------------------------
inline void FloatCompressor::GetMinMax(F32& pMin, F32& pMax)
{
  fprintf(stderr,"not yet implemented\n");
}
//-----------------------------------------------------------------------------
#endif
