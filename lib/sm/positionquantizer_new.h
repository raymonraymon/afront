/*
===============================================================================

  FILE:  positionquantizer_new.h
  
  CONTENTS:
  
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    29 July 2004 -- adapted from the old PositionQuantizer. this one is better
  
===============================================================================
*/
#ifndef POSITION_QUANTIZER_NEW_H
#define POSITION_QUANTIZER_NEW_H

#include "mydefs.h"

class PositionQuantizerNew
{
public:
  // Reset:
  inline void Reset();

  // SetPrecision:
  inline void SetPrecision(U32 uBits);
  // GetPrecision:
  inline void GetPrecision(U32 &uBits);

  // SetRange:
  inline void SetRange(U32 uRange);
  // GetRange:
  inline void GetRange(U32 &uRange);

  // UpdateMinMax:
  inline void UpdateMinMax(F32* pPos);
  // SetMinMax:
  inline void SetMinMax(F32* pMin, F32* pMax);
  // GetMinMax:
  inline void GetMinMax(F32* pMin, F32* pMax);

  // SetupQuantizer:
  inline void SetupQuantizer();

  // EnQuantize:
  inline void EnQuantize(const F32* pPosition, I32* pQuantPosition);
  // DeQuantize:
  inline void DeQuantize(const I32* pQuantPosition, F32* pPosition);

  // Clamp:
  inline void Clamp(I32* pQuantPosition);
  inline void Clamp(F32* pPosition);
  inline void Clamp(const F32* pPosition, I32* pQuantPosition);

  // Constructor:
  PositionQuantizerNew();
  // Destructor:
  ~PositionQuantizerNew();

//private:
  // Clamp:
  inline I32 Clamp(I32 iValue, I32 iMin, I32 iMax);
  inline F32 Clamp(F32 fValue, F32 fMin, F32 fMax);
  inline I32 Clamp(F32 fValue, I32 iMin, I32 iMax);
  
  U32 m_uBits;
  U32 m_uRange;

  I32 m_aiQuantMin[3];
  I32 m_aiQuantMax[3];
  I32 m_aiQuantRange[3];
  
  F32 m_afMin[3];
  F32 m_afMax[3];

  F64 m_dEnQuantizeMultiplier;
  F64 m_dDeQuantizeMultiplier;

}; // end class PositionQuantizerNew

//-----------------------------------------------------------------------------
//  Public Inline methods
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Reset:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::Reset()
{
  m_uBits = 12;
  m_uRange = 1 << 12;
  m_afMin[0] = m_afMin[1] = m_afMin[2] = F32_MAX;
  m_afMax[0] = m_afMax[1] = m_afMax[2] = F32_MIN;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SetPrecision:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::SetPrecision(U32 uBits)
{
  m_uBits = uBits;
  m_uRange = 1 << uBits;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GetPrecision:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::GetPrecision(U32 &uBits)
{
  uBits = m_uBits;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SetPrecision:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::SetRange(U32 uRange)
{
  m_uRange = uRange;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GetRange:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::GetRange(U32 &uRange)
{
  uRange = m_uRange;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// UpdateMinMax:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::UpdateMinMax(F32* pPos)
{
  if (pPos != NULL)
  {
    if (pPos[0] < m_afMin[0]) m_afMin[0] = pPos[0];
    if (pPos[1] < m_afMin[1]) m_afMin[1] = pPos[1];
    if (pPos[2] < m_afMin[2]) m_afMin[2] = pPos[2];
    if (pPos[0] > m_afMax[0]) m_afMax[0] = pPos[0];
    if (pPos[1] > m_afMax[1]) m_afMax[1] = pPos[1];
    if (pPos[2] > m_afMax[2]) m_afMax[2] = pPos[2];
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SetMinMax:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::SetMinMax(F32* pMin, F32* pMax)
{
  if ((pMin != NULL) && (pMax != NULL))
  {
    m_afMin[0] = pMin[0];
    m_afMin[1] = pMin[1];
    m_afMin[2] = pMin[2];
    m_afMax[0] = pMax[0];
    m_afMax[1] = pMax[1];
    m_afMax[2] = pMax[2];
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GetMinMax:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::GetMinMax(F32* pMin, F32* pMax)
{
  if ((pMin != NULL) && (pMax != NULL))
  {
    pMin[0] = m_afMin[0];
    pMin[1] = m_afMin[1];
    pMin[2] = m_afMin[2];
    pMax[0] = m_afMax[0];
    pMax[1] = m_afMax[1];
    pMax[2] = m_afMax[2];
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SetupQuantizer:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::SetupQuantizer()
{
  I32 iNumSteps = m_uRange - 1;
  F64 dMaxRange;
  F64 dRange0 = ((F64)m_afMax[0]-(F64)m_afMin[0]); // is that casting a good idea? jack says ... it does not do anything
  F64 dRange1 = ((F64)m_afMax[1]-(F64)m_afMin[1]);
  F64 dRange2 = ((F64)m_afMax[2]-(F64)m_afMin[2]);
  if (dRange0 > dRange1)
  {
    if (dRange0 > dRange2)
    {
      dMaxRange = dRange0;
    }
    else
    {
      dMaxRange = dRange2;
    }
  }
  else
  {
    if (dRange1 > dRange2)
    {
      dMaxRange = dRange1;
    }
    else
    {
      dMaxRange = dRange2;
    }
  }

  if (dMaxRange > 0.0)
  {
    m_dEnQuantizeMultiplier = (F64) iNumSteps / dMaxRange; // is that casting a good idea? jack says ... it's automatically promoted
    m_dDeQuantizeMultiplier = dMaxRange / (F64) iNumSteps;
    m_aiQuantMin[0] = 0;
    m_aiQuantMin[1] = 0;
    m_aiQuantMin[2] = 0;
    m_aiQuantMax[0] = I32(dRange0*m_dEnQuantizeMultiplier + 0.5);
    m_aiQuantMax[1] = I32(dRange1*m_dEnQuantizeMultiplier + 0.5);
    m_aiQuantMax[2] = I32(dRange2*m_dEnQuantizeMultiplier + 0.5);
    m_aiQuantRange[0] = m_aiQuantMax[0] - m_aiQuantMin[0] + 1;
    m_aiQuantRange[1] = m_aiQuantMax[1] - m_aiQuantMin[1] + 1;
    m_aiQuantRange[2] = m_aiQuantMax[2] - m_aiQuantMin[2] + 1;
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// EnQuantize:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::EnQuantize(const F32* pPosition, I32* pQuantPosition)
{
  pQuantPosition[0] = I32(m_dEnQuantizeMultiplier * ((F64)pPosition[0] - (F64)m_afMin[0]) + 0.5); // is that casting a good idea? jack says ... it is important to avoid bit cancellation
  pQuantPosition[1] = I32(m_dEnQuantizeMultiplier * ((F64)pPosition[1] - (F64)m_afMin[1]) + 0.5);
  pQuantPosition[2] = I32(m_dEnQuantizeMultiplier * ((F64)pPosition[2] - (F64)m_afMin[2]) + 0.5);

  Clamp(pQuantPosition); // just to be sure
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// DeQuantize:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::DeQuantize(const I32* pQuantPosition, F32* pPosition)
{
  pPosition[0] = (F32)(m_dDeQuantizeMultiplier*pQuantPosition[0] + (F64)m_afMin[0]);// is that casting a good idea? jack says ... yes.
  pPosition[1] = (F32)(m_dDeQuantizeMultiplier*pQuantPosition[1] + (F64)m_afMin[1]);
  pPosition[2] = (F32)(m_dDeQuantizeMultiplier*pQuantPosition[2] + (F64)m_afMin[2]);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Clamp:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::Clamp(I32* pQuantPosition)
{
  pQuantPosition[0] = Clamp(pQuantPosition[0], m_aiQuantMin[0], m_aiQuantMax[0]);
  pQuantPosition[1] = Clamp(pQuantPosition[1], m_aiQuantMin[1], m_aiQuantMax[1]);
  pQuantPosition[2] = Clamp(pQuantPosition[2], m_aiQuantMin[2], m_aiQuantMax[2]);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Clamp:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::Clamp(F32* pPosition)
{
  pPosition[0] = Clamp(pPosition[0], m_afMin[0], m_afMax[0]);
  pPosition[1] = Clamp(pPosition[1], m_afMin[1], m_afMax[1]);
  pPosition[2] = Clamp(pPosition[2], m_afMin[2], m_afMax[2]);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Clamp:
//-----------------------------------------------------------------------------
inline void PositionQuantizerNew::Clamp(const F32* pPosition, I32* pQuantPosition)
{
  pQuantPosition[0] = Clamp(pPosition[0], m_aiQuantMin[0], m_aiQuantMax[0]);
  pQuantPosition[1] = Clamp(pPosition[1], m_aiQuantMin[1], m_aiQuantMax[1]);
  pQuantPosition[2] = Clamp(pPosition[2], m_aiQuantMin[2], m_aiQuantMax[2]);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor:
//-----------------------------------------------------------------------------
inline PositionQuantizerNew::PositionQuantizerNew()
{
  Reset();
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Destructor:
//-----------------------------------------------------------------------------
inline PositionQuantizerNew::~PositionQuantizerNew()
{
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//  Private Inline methods
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Clamp:
//-----------------------------------------------------------------------------
inline I32 PositionQuantizerNew::Clamp(I32 iValue, I32 iMin, I32 iMax)
{
  return (iValue<=iMin) ? iMin : (iValue>=iMax) ? iMax : iValue;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Clamp:
//-----------------------------------------------------------------------------
inline F32 PositionQuantizerNew::Clamp(F32 fValue, F32 fMin, F32 fMax)
{
  return (fValue<=fMin) ? fMin : (fValue>=fMax) ? fMax : fValue;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Clamp:
//-----------------------------------------------------------------------------
inline I32 PositionQuantizerNew::Clamp(F32 fValue, I32 iMin, I32 iMax)
{
  return (fValue<=iMin) ? iMin : (fValue>=iMax) ? iMax : (I32)(fValue+0.5f);
}
//-----------------------------------------------------------------------------

#endif
