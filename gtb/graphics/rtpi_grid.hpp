
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


#ifndef GTB_RTPI_GRID_INCLUDED
#define GTB_RTPI_GRID_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/model_rtpi.hpp>


GTB_BEGIN_NAMESPACE


class RtpiGrid {
public:
	explicit RtpiGrid(const ModelRtpi &model);
  	~RtpiGrid();

	unsigned num_rows() const;
	unsigned num_cols() const;
	const Rtpi &sample(unsigned i, unsigned j) const;

protected:

	Rtpi **_data;
	unsigned _num_rows;
	unsigned _num_cols;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/rtpi_grid.ipp>
#endif

#endif // GTB_RTPI_GRID_INCLUDED
