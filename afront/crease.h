
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


#ifndef __CREASE_H_CSCHEID_20051213
#define __CREASE_H_CSCHEID_20051213

#include <vector>
#include <map>
#include <list>
#include "common.h"
#include "parallel.h"
#include "guidance.h"
#include "triangulate_mesh.h"

class CreaseExtractor
{
    // Terminology: a "kink" is formed by a set of "creases". Each
    // kink is a subset of all the crease edges between two kink
    // vertices. The kink vertices are vertices with crease degree >2,
    // and every vertex in the middle of the kink has degree 2.
    // Exceptionally, there might exist loop kinks, which are a loop
    // of crease edges all with degree two (think cylinder edges). We
    // arbitrarily choose a vertex of the loop to be a kink vertex,
    // and the kink graph will have a self-edge.  There might also
    // exist "twig" kinks (twig as in a short branch), which are kinks
    // that vanish before connecting to some other kink vertex. These
    // are handled by inserting a kink vertex artificially at the end
    // of the twig kink, thus splitting the degenerate kink into two
    // regular kinks.
    
public:
    CreaseExtractor(const TriangleMesh &_m,
		    MeshGuidanceField *_guidance,
		    MeshProjector *_projector,
		    std::vector<gtb::Point3> &cp,
		    std::vector<gtb::Vector3> &con,
		    std::vector< std::vector<int> > &ci,
		    std::vector< std::vector<gtb::Vector3> > &cn,
		    std::vector<int> &corner_tris,
		    float s = -2.0,
		    unsigned _small_crease_threshold = 0):
	m(_m),
	guidance(_guidance),
	projector(_projector),
	crease_points(cp),
	crease_onormals(con),
	crease_indices(ci),
	crease_normals(cn),
	corner_triangles(corner_tris),
	sharp(s),
	small_crease_threshold(_small_crease_threshold)
    {};

    void Extract();

private:
    typedef std::map<int, std::list<int> > Graph;

    struct ResampledCurve {
	std::vector<gtb::Point3> points;
	std::vector<gtb::Vector3> normals;
	std::vector<int> indices;
    };

    const TriangleMesh &m;
    MeshGuidanceField *guidance;
    MeshProjector *projector;
    std::vector<gtb::Point3> &crease_points;
    std::vector<gtb::Vector3> &crease_onormals;
    std::vector< std::vector<int> > &crease_indices;
    std::vector< std::vector<gtb::Vector3> > &crease_normals;
    std::vector<int> &corner_triangles;
    float sharp;
    unsigned small_crease_threshold;

    Graph creases;
    std::list<int> kink_vertices;
    std::vector<std::vector<int> > kinks;
    std::map< std::pair<int, int>, ResampledCurve > resampledKinks;
    Graph kink_graph;
    int total_half_edges;

    void MakeGraph();
    bool IsKinkVertex(int v);
    int FindKinkVertex(const Graph &c);
    void ComputeKinks();
    gtb::Vector3 ComputeEdgeNormal(std::vector<int>&, int);
    void ComputeResampledKinks();
    int FindKinkWithEndPoint(int v, int going_through=-1);
    void AddCornerTriangles();

    // Illustration of OtherKinkEndPoint and KinkVertexNeighbor

    //  /- endpoint
    // v    v-> KinkVertexNeighbor(kink_index, endpoint)       v- OtherKinkEndPoint(kink_index, endpoint)
    // *----*----*----*----*-- ... -- *----*----*----*----*----*
    // |                                                       |
    // \-------------------------------------------------------/
    //            The kink (kinks[kink_index])

    void DeleteSmallCreases();
    int OtherKinkEndPoint(int kink_index, int endpoint);
    int KinkVertexNeighbor(int kink_index, int endpoint);
    void AddKinkToFront(std::vector<int> &front, std::vector<gtb::Vector3> &normals, int kink_index, int endpoint);
    void ComputeFronts();
    void SetOutput();
    void ReverseFrontIfNecessary(std::vector<int> &front, std::vector<gtb::Vector3> &normals);
};

#endif
