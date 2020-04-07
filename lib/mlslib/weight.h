
/*

Copyright 2007 University of Utah


This file is part of Afront.

Afront is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Afront is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA

*/


#ifndef __CORNER_WEIGHTS_H
#define __CORNER_WEIGHTS_H

#include "mlslibdefs.h"

MLSLIB_BEGIN_NAMESPACE

template <typename REAL>
class WeightFunction
{
public:
    typedef REAL areal;

    WeightFunction() {}
    virtual ~WeightFunction() {}

    //
    // Compute the weight for the input value
    // The input is value^2, since it is commonly the
    // distance and the squared value is probably what we
    // need, thus there is no reason to compute the square of
    // the square root.
    virtual areal weight(areal squared_value) const = 0;

    //
    // Set the "scale" of the function
    // by default the function scale is 1
    // i.e. it is defined over the range [0,1]
    // (or have some other notion of scale with value of 1)
    //
    virtual void set_radius(areal radius) = 0;

    //
    // This function returns the "scale" of the function
    // If the function is defined over [0, infinity), then
    // this function returns the scale for the given percentile
    // i.e the \int_0^{percentile} / \int_0^\inf >= percentile
    //
    virtual areal influence_radius(areal percentile) const = 0;

    //
    // cloning
    //
    virtual WeightFunction* clone() const = 0;
};

template <typename REAL>
class GaussianWeightFunction : public WeightFunction<REAL>
{
public:
    GaussianWeightFunction() {}
    GaussianWeightFunction(typename WeightFunction<REAL>::areal sigma) { set_radius(sigma); }
    GaussianWeightFunction(const GaussianWeightFunction& rhs) : _sigma(rhs._sigma), _factor(rhs._factor) {}

    typename WeightFunction<REAL>::areal 
    weight(typename WeightFunction<REAL>::areal squared_value) const 
    {
        return exp(squared_value * _factor); 
    }

    void set_radius(typename WeightFunction<REAL>::areal sigma)
    { 
        _sigma = sigma; 
        _factor = -1.0 / (2*sigma*sigma);
    }

    typename WeightFunction<REAL>::areal 
    influence_radius(typename WeightFunction<REAL>::areal percentile) const
    {
// BUGBUG, this is not what we mean
        return 2*_sigma;
    }

    WeightFunction<REAL>* clone() const { return new GaussianWeightFunction(*this); }

    // protected:
    typename WeightFunction<REAL>::areal _sigma;
    typename WeightFunction<REAL>::areal _factor;
};


/*--------------------- [ bounded weight 4 ] ------------------------------------*/
//
// The following function that has only powers of 2, so it's easy
// to compute the derivative of for sqrt(x)
// 1 - 2x^2 + x^4
//

template <typename REAL>
class BoundedWeightFunction4 : public WeightFunction<REAL>
{
public:
    BoundedWeightFunction4() {}
    BoundedWeightFunction4(typename WeightFunction<REAL>::areal radius) { set_radius(radius); }
    BoundedWeightFunction4(const BoundedWeightFunction4& rhs) : _radius(rhs._radius), _radius2(rhs._radius2), _factor(rhs._factor) {}

    typename WeightFunction<REAL>::areal weight(typename WeightFunction<REAL>::areal squared_value) const;

    void set_radius(typename WeightFunction<REAL>::areal radius)
    { 
        _radius = radius;
        _radius2 = _radius*_radius;
        _factor = 1.0 / _radius2;
    }

    typename WeightFunction<REAL>::areal influence_radius(typename WeightFunction<REAL>::areal percentile) const
    {
        return _radius;
    }

    WeightFunction<REAL>* clone() const { return new BoundedWeightFunction4(*this); }

    // protected:
    typename WeightFunction<REAL>::areal _radius;
    typename WeightFunction<REAL>::areal _radius2;
    typename WeightFunction<REAL>::areal _factor;
protected:
};

template <typename REAL>
inline 
typename WeightFunction<REAL>::areal
BoundedWeightFunction4<REAL>::weight(typename WeightFunction<REAL>::areal squared_value) const 
{
    if (squared_value > _radius2) return 0;
    else 
    {
        squared_value *= _factor;
        return 1 - ((typename WeightFunction<REAL>::areal)2)*squared_value + squared_value*squared_value;
    }
}

/*--------------------- [ bounded weight 4 ] ------------------------------------*/

template <typename REAL>
class BoundedWeightFunction : public WeightFunction<REAL>
{
public:
    BoundedWeightFunction() {}
    BoundedWeightFunction(typename WeightFunction<REAL>::areal radius) { set_radius(radius); }
    BoundedWeightFunction(const BoundedWeightFunction& rhs) : _radius(rhs._radius), _factor(rhs._factor) {}

    typename WeightFunction<REAL>::areal weight(typename WeightFunction<REAL>::areal squared_value) const;

    void set_radius(typename WeightFunction<REAL>::areal radius)
    { 
        _radius = radius;
        _factor = 1.0 / radius;
    }

    typename WeightFunction<REAL>::areal influence_radius(typename WeightFunction<REAL>::areal percentile) const
    {
        return _radius;
    }

    WeightFunction<REAL>* clone() const { return new BoundedWeightFunction(*this); }

    // protected:
    typename WeightFunction<REAL>::areal _radius;
    typename WeightFunction<REAL>::areal _factor;
protected:
    REAL lsqrt(REAL val) const;
};

template<> inline float BoundedWeightFunction<float>::lsqrt(float val) const {return sqrtf(val); }
template<> inline double BoundedWeightFunction<double>::lsqrt(double val) const {return sqrt(val); }

template <typename REAL>
inline 
typename WeightFunction<REAL>::areal
BoundedWeightFunction<REAL>::weight(typename WeightFunction<REAL>::areal squared_value) const 
{
//    return gtb::bounded_decay(lsqrt(squared_value) * _factor);
    return gtb::infinit_bounded_decay(lsqrt(squared_value) * _factor);
}

MLSLIB_END_NAMESPACE

#endif // __CORNER_WEIGHTS_H
