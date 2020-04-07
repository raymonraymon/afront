
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
#ifndef WIN32
#include <gtb/graphics/model_rtpi.hpp>
#include <gtb/io/io.hpp>
#include <gtb/math/math.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/model_rtpi.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


ModelRtpi::ModelRtpi()
{
}


ModelRtpi::~ModelRtpi()
{
	for (unsigned i = 0; i < _strips.size(); i++) {
		delete _strips[i];
		_strips[i] = NULL;
	}
}


void ModelRtpi::read_header(FILE *fp,
			    unsigned &n_strips,
			    unsigned &a_max_points)
{
	assert(NULL != fp);
	char header[80];
    do
    {
	    fgets(header, sizeof(header), fp);
    } while ((!feof(fp)) && (header[0] == '#'));
	unsigned dummy;
	if (sscanf(header, "P9 %u %u", &dummy, &dummy) == 2) {
		GTB_ERROR("RTPI version 1 not supported");
	} else if (sscanf(header, "RTPIv2 %u %u",
			  &n_strips, &a_max_points) != 2) {
		GTB_ERROR("unrecognized file format");
	}
}


void ModelRtpi::read(FILE *fp)
{
	unsigned n_strips;
	read_header(fp, n_strips, _max_points);

	_min_r = real::MAX;
	_max_r = -real::MAX;

	_min_t = 2 * M_PI;
	_max_t = -2 * M_PI;

	_min_p = 2 * M_PI;
	_max_p = -2 * M_PI;

	_min_i = INT_MAX;
	_max_i = INT_MIN;

	for (unsigned i = 0; i < n_strips; i++) {
		_strips.push_back(new RtpiStrip);
		unsigned n;
		read_unsigned(&n, fp);
		for (unsigned j = 0; j < n; j++) {
			strip(i).add_point(Rtpi());
			strip(i)[j].read(fp);

			_min_r = min(_min_r, strip(i)[j].r());
			_max_r = max(_max_r, strip(i)[j].r());

			_min_t = min(_min_t, strip(i)[j].t());
			_max_t = max(_max_t, strip(i)[j].t());

			_min_p = min(_min_p, strip(i)[j].p());
			_max_p = max(_max_p, strip(i)[j].p());

			_min_i = min(_min_i, strip(i)[j].i());
			_max_i = max(_max_i, strip(i)[j].i());
		}
	}
}


void ModelRtpi::write(FILE *fp)
{
	assert(NULL != fp);

	fprintf(fp, "RTPIv2     %d  %d\n", num_strips(), max_points());
	for (unsigned i = 0; i < num_strips(); i++) {
		write_unsigned(strip(i).num_points(), fp);
		for (unsigned j = 0; j < strip(i).num_points(); j++) {
			strip(i)[j].write(fp);
		}
	}
}


void ModelRtpi::compute_bounding_box() const
{
	// should this be optimized?
	vector<Point3> points;
	for (unsigned i = 0; i < num_strips(); i++) {
		for (unsigned j = 0; j < strip(i).num_points(); j++) {
			points.push_back(strip(i)[j].point());
		}
	}
	_bounding_box = Box3::bounding_box(points);
}


void ModelRtpi::compute_centroid() const
{
	// should this be optimized?
	vector<Point3> points;
	for (unsigned i = 0; i < num_strips(); i++) {
		for (unsigned j = 0; j < strip(i).num_points(); j++) {
			points.push_back(strip(i)[j].point());
		}
	}
	_centroid = Point3::centroid(points);
}


void ModelRtpi::histogram(vector<unsigned> &hist)
{
	// init histogram
	for (unsigned i = 0; i < hist.size(); i++) {
		hist[i] = 0;
	}

	histogram_add(hist);
}


void ModelRtpi::histogram_add(vector<unsigned> &hist)
{
	// for each strip
	for (unsigned i = 0; i < num_strips(); i++) {
		const RtpiStrip &s = strip(i);

		// for each point
		for (unsigned j = 0; j < s.num_points(); j++) {
			const Rtpi &p = s[j];

			// update histogram
  			real_t c = (real_t) p.i() / (real_t) max_intensity();
			hist[(unsigned) (c * (hist.size() - 1))]++;
		}
	}
}


void ModelRtpi::equalize_histogram(unsigned hist_size)
{
	assert(hist_size > 0);

	// compute histogram
	vector<unsigned> hist(hist_size);
	histogram(hist);

	// compute cumulative probabilities
	unsigned n = num_points();
	vector<real_t> cumulative_prob(hist_size);
	cumulative_prob[0] = (real_t) hist[0] / (real_t) n;
	for (unsigned i = 1; i < hist_size; i++) {
		cumulative_prob[i] = (cumulative_prob[i - 1]
				      + (real_t) hist[i] / (real_t) n);
	}

	// update intensities
	for (unsigned i = 0; i < num_strips(); i++) {
		RtpiStrip &s = strip(i);
		for (unsigned j = 0; j < s.num_points(); j++) {
			Rtpi &p = s[j];
			real_t c = (real_t) p.i() / (real_t) max_intensity();
			unsigned k = (unsigned) (c * (hist.size() - 1));
			p.set_i((int) (cumulative_prob[k] * max_intensity()));
		}
	}
}


void ModelRtpi::scale(real_t a)
{
	for (unsigned i = 0; i < num_strips(); i++) {
		RtpiStrip &s = strip(i);
		for (unsigned j = 0; j < s.num_points(); j++) {
			Rtpi &p = s[j];
			p.set_r(p.r() * a);
		}
	}
}

void ModelRtpi::insert_strip(const RtpiStrip& a_strip)
{
	_strips.push_back(new RtpiStrip(a_strip));
}

GTB_END_NAMESPACE
