
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


#include <gtb/gtb.hpp>
#include "surfelhierarchy.h"
#ifndef WIN32
#include "surfel_set_io.hpp"
#include <gtb/io/io.hpp>
#include <gtb/graphics/ogltools.h>
#endif // WIN32

GTB_BEGIN_NAMESPACE

void surfel_hierarchy::render() const
{
    render(0, _hierarchy.size());
}

void surfel_hierarchy::render(unsigned a_level) const
{
    assert(a_level < _hierarchy.size());
    const surfel_set& surfelset = *_hierarchy[a_level];
    unsigned n_points = surfelset.size();
    glBegin(GL_POINTS);
    for (unsigned i = 0; i < n_points; ++i)
    {
        surfelset.normal(i).load_as_normal();
        surfelset.vertex(i).load();
    }
    glEnd();
}

void surfel_hierarchy::render(unsigned level1, unsigned level2) const
{
    for (unsigned i = level1; i < level2; ++i)
    {
        render(i);
    }
}

void surfel_hierarchy::write_as_obj(const char* name) const
{
    afree<FILE*> f(xfopen(name, "w"), fclose);

    for (unsigned i = 0; i < num_levels(); ++i)
    {
		write_obj(f, *(_hierarchy[i]));
    }
}

unsigned surfel_hierarchy::size_() const
{
	unsigned s = 0;
    for (unsigned i = 0; i < num_levels(); ++i)
    {
		s += _hierarchy[i]->size();
    }
	return s;
}

void surfel_hierarchy::clear_level(unsigned level_)
{
	assert(level_ < num_levels());
	_hierarchy[level_]->clear();
}

GTB_END_NAMESPACE


