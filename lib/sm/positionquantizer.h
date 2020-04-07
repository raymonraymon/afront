/*
===============================================================================

  FILE:  positionquantizer.h
  
  CONTENTS:
  
    Triangle Strip Compression - GI 2000 demonstration software
  
  PROGRAMMERS:
  
    Martin Isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    Copyright (C) 2000  Martin Isenburg (isenburg@cs.unc.edu)
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    28 June 2000 -- Martin - finalized for GI submission.
  
===============================================================================
*/
#ifndef POSITION_QUANTIZER_H
#define POSITION_QUANTIZER_H

#include "mydefs.h"

class PositionQuantizer
{
public:
  // Reset:
  inline void Reset();

  // SetPrecision:
  inline void SetPrecision(U32 uBits);
  // GetPrecision:
  inline void GetPrecision(U32 &uBits);

  // UpdateMinMax:
  inline void UpdateMinMax(F32* pPos);
  // SetMinMax:
  inline void SetMinMax(F32* pMin, F32* pMax);
  // GetMinMax:
  inline void GetMinMax(F32* pMin, F32* pMax);

  // SetupQuantizer:
  inline void SetupQuantizer();

  // EnQuantize:
  inline void EnQuantize(F32* pPosition, I32* pQuantPosition);
  // DeQuantize:
  inline void DeQuantize(I32* pQuantPosition, F32* pPosition);

  // Clamp:
  inline void Clamp(I32* pQuantPosition);
  // Wrap:
  inline void Wrap(I32* pQuantPosition);

  // GetCorrector:
  inline void GetCorrector(I32* pPredPosition, I32* pRealPosition, I32* pCorrector);
  // AddCorrector:
  inline void AddCorrector(I32* pPredPosition, I32* pCorrector, I32* pRealPosition);

  // Constructor:
  PositionQuantizer();
  // Destructor:
  ~PositionQuantizer();

//private:
  // Clamp:
  inline I32 Clamp(I32 iValue, I32 iMin, I32 iMax);
  // Wrap:
  inline I32 Wrap(I32 iValue, I32 iMin, I32 iMax, I32 iRange);
  
  U32 m_uBits;

  I32 m_aiMinCode[3];
  I32 m_aiMaxCode[3];

  I32 m_aiRangeCode[3];

  I32 m_aiMinCorrector[3];
  I32 m_aiMaxCorrector[3];

  I32 m_aiAbsRangeCorrector[3];
  
  F32 m_afMin[3];
  F32 m_afMax[3];

  F64 m_dEnQuantizeMultiplier;
  F64 m_dDeQuantizeMultiplier;

}; // end class PositionQuantizer

//-----------------------------------------------------------------------------
//  Public Inline methods
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Reset:
//-----------------------------------------------------------------------------
inline void PositionQuantizer::Reset()
{
  m_uBits = 12;
  m_afMin[0] = m_afMin[1] = m_afMin[2] = F32_MAX;
  m_afMax[0] = m_afMax[1] = m_afMax[2] = F32_MIN;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SetPrecision:
//-----------------------------------------------------------------------------
inline void PositionQuantizer::SetPrecision(U32 uBits)
{
  if ((0 < uBits) && (uBits <= 30))
  {
    m_uBits = uBits;
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GetPrecision:
//-----------------------------------------------------------------------------
inline void PositionQuantizer::GetPrecision(U32 &uBits)
{
  uBits = m_uBits;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// UpdateMinMax:
//-----------------------------------------------------------------------------
inline void PositionQuantizer::UpdateMinMax(F32* pPos)
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
inline void PositionQuantizer::SetMinMax(F32* pMin, F32* pMax)
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
inline void PositionQuantizer::GetMinMax(F32* pMin, F32* pMax)
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
inline void PositionQuantizer::SetupQuantizer()
{
  I32 iNumSteps = (1 << m_uBits) - 2; // will be changed to -1 soon
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
    m_aiMinCode[0] = 0;
    m_aiMinCode[1] = 0;
    m_aiMinCode[2] = 0;
    m_aiMaxCode[0] = I32(dRange0*m_dEnQuantizeMultiplier + 0.5);
    m_aiMaxCode[1] = I32(dRange1*m_dEnQuantizeMultiplier + 0.5);
    m_aiMaxCode[2] = I32(dRange2*m_dEnQuantizeMultiplier + 0.5);
    m_aiRangeCode[0] = m_aiMaxCode[0] - m_aiMinCode[0] + 1;
    m_aiRangeCode[1] = m_aiMaxCode[1] - m_aiMinCode[1] + 1;
    m_aiRangeCode[2] = m_aiMaxCode[2] - m_aiMinCode[2] + 1;
    m_aiMinCorrector[0] = m_aiMinCode[0]-(m_aiRangeCode[0]/2);
    m_aiMinCorrector[1] = m_aiMinCode[1]-(m_aiRangeCode[1]/2);
    m_aiMinCorrector[2] = m_aiMinCode[2]-(m_aiRangeCode[2]/2);
    m_aiMaxCorrector[0] = m_aiMaxCode[0]-(m_aiRangeCode[0]/2);
    m_aiMaxCorrector[1] = m_aiMaxCode[1]-(m_aiRangeCode[1]/2);
    m_aiMaxCorrector[2] = m_aiMaxCode[2]-(m_aiRangeCode[2]/2);
    m_aiAbsRangeCorrector[0] = (m_aiRangeCode[0])/2+1;
    m_aiAbsRangeCorrector[1] = (m_aiRangeCode[1])/2+1;
    m_aiAbsRangeCorrector[2] = (m_aiRangeCode[2])/2+1;
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// EnQuantize:
//-----------------------------------------------------------------------------
inline void PositionQuantizer::EnQuantize(F32* pPosition, I32* pQuantPosition)
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
inline void PositionQuantizer::DeQuantize(I32* pQuantPosition, F32* pPosition)
{
  pPosition[0] = (F32)(m_dDeQuantizeMultiplier*pQuantPosition[0] + (F64)m_afMin[0]);// is that casting a good idea? jack says ... yes.
  pPosition[1] = (F32)(m_dDeQuantizeMultiplier*pQuantPosition[1] + (F64)m_afMin[1]);
  pPosition[2] = (F32)(m_dDeQuantizeMultiplier*pQuantPosition[2] + (F64)m_afMin[2]);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Clamp:
//-----------------------------------------------------------------------------
inline void PositionQuantizer::Clamp(I32* pQuantPosition)
{
  pQuantPosition[0] = Clamp(pQuantPosition[0], m_aiMinCode[0], m_aiMaxCode[0]);
  pQuantPosition[1] = Clamp(pQuantPosition[1], m_aiMinCode[1], m_aiMaxCode[1]);
  pQuantPosition[2] = Clamp(pQuantPosition[2], m_aiMinCode[2], m_aiMaxCode[2]);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GetCorrector:
//-----------------------------------------------------------------------------
inline void PositionQuantizer::GetCorrector(I32* pPredPosition, I32* pRealPosition, I32* pCorrector)
{
  if (pPredPosition == NULL) // then use default predictor [0/0/0]
  {
    pCorrector[0] = Wrap(pRealPosition[0], m_aiMinCorrector[0], m_aiMaxCorrector[0], m_aiRangeCode[0]);
    pCorrector[1] = Wrap(pRealPosition[1], m_aiMinCorrector[1], m_aiMaxCorrector[1], m_aiRangeCode[1]);
    pCorrector[2] = Wrap(pRealPosition[2], m_aiMinCorrector[2], m_aiMaxCorrector[2], m_aiRangeCode[2]);
  }
  else
  {
    pCorrector[0] = Wrap(pRealPosition[0] - pPredPosition[0], m_aiMinCorrector[0], m_aiMaxCorrector[0], m_aiRangeCode[0]);
    pCorrector[1] = Wrap(pRealPosition[1] - pPredPosition[1], m_aiMinCorrector[1], m_aiMaxCorrector[1], m_aiRangeCode[1]);
    pCorrector[2] = Wrap(pRealPosition[2] - pPredPosition[2], m_aiMinCorrector[2], m_aiMaxCorrector[2], m_aiRangeCode[2]);
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// AddCorrector:
//-----------------------------------------------------------------------------
inline void PositionQuantizer::AddCorrector(I32* pPredPosition, I32* pCorrector, I32* pRealPosition)
{
  if (pPredPosition == NULL) // then use default predictor [0/0/0]
  {
    pRealPosition[0] = Wrap(pCorrector[0], m_aiMinCode[0], m_aiMaxCode[0], m_aiRangeCode[0]);
    pRealPosition[1] = Wrap(pCorrector[1], m_aiMinCode[1], m_aiMaxCode[1], m_aiRangeCode[1]);
    pRealPosition[2] = Wrap(pCorrector[2], m_aiMinCode[2], m_aiMaxCode[2], m_aiRangeCode[2]);
  }
  else
  {
    pRealPosition[0] = Wrap(pPredPosition[0] + pCorrector[0], m_aiMinCode[0], m_aiMaxCode[0], m_aiRangeCode[0]);
    pRealPosition[1] = Wrap(pPredPosition[1] + pCorrector[1], m_aiMinCode[1], m_aiMaxCode[1], m_aiRangeCode[1]);
    pRealPosition[2] = Wrap(pPredPosition[2] + pCorrector[2], m_aiMinCode[2], m_aiMaxCode[2], m_aiRangeCode[2]);
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor:
//-----------------------------------------------------------------------------
inline PositionQuantizer::PositionQuantizer()
{
  Reset();
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Destructor:
//-----------------------------------------------------------------------------
inline PositionQuantizer::~PositionQuantizer()
{
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//  Private Inline methods
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Clamp:
//-----------------------------------------------------------------------------
inline I32 PositionQuantizer::Clamp(I32 iValue, I32 iMin, I32 iMax)
{
  return (iValue<=iMin) ? iMin : (iValue>=iMax) ? iMax : iValue;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Wrap:
//-----------------------------------------------------------------------------
inline I32 PositionQuantizer::Wrap(I32 iValue, I32 iMin, I32 iMax, I32 iRange)
{
  while (iValue<iMin)
  {
    iValue += iRange;
  }
  while (iValue>iMax)
  {
    iValue -= iRange;
  }
  return iValue;
}
//-----------------------------------------------------------------------------

#endif
