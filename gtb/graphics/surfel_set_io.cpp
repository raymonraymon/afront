
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
#include <gtb/memory/ptrs.h>

#include "surfel_set_io.hpp"

#ifndef WIN32
#include <gtb/io/io.hpp>
#endif // WIN32

GTB_BEGIN_NAMESPACE

template<class T>
void read_points(const char* name, tsurfel_set<T>& surfels)
{
    const char* pdot = strrchr(name, '.');
    if (pdot == 0) return;
    if (stricmp(pdot+1, "obj") == 0) read_obj(name, surfels);
    else if (stricmp(pdot+1, "sst") == 0) read_bin_points(name, surfels);
    else return;
}

template void read_points(const char* name, tsurfel_set<float>& surfels);
template void read_points(const char* name, tsurfel_set<double>& surfels);


template<class T>
void write_points(const char* name, const tsurfel_set<T>& surfels)
{
    const char* pdot = strrchr(name, '.');
    if (pdot == 0) return;
    if (stricmp(pdot+1, "obj") == 0) write_obj(name, surfels);
    else if (stricmp(pdot+1, "sst") == 0) write_bin_points(name, surfels);
    else return;
}

template void write_points(const char* name, const tsurfel_set<float>& surfels);
template void write_points(const char* name, const tsurfel_set<double>& surfels);

template<class T>
void append_points(const char* name, const tsurfel_set<T>& surfels)
{
    const char* pdot = strrchr(name, '.');
    if (pdot == 0) return;
    if (stricmp(pdot+1, "obj") == 0)
    {
        afree<FILE*> f(xfopen(name, "a"), fclose);
        fseek(f, 0, SEEK_END);
        write_obj(f, surfels);
    }
    else if (stricmp(pdot+1, "sst") == 0) append_bin_points(name, surfels);
    else return;
}

template void append_points(const char* name, const tsurfel_set<float>& surfels);
template void append_points(const char* name, const tsurfel_set<double>& surfels);

template<class T>
void read_obj(const char* name, tsurfel_set<T>& surfels)
{
    afree<FILE*> f(xfopen(name, "r"), fclose);

//    surfels.reserve(6600000, true, true, false);

    while (!feof(f))
    {
        char line[1000];
        if (fgets(line, sizeof(line), f) == 0) break;


        if (strncmp(line, "v ", 2) == 0)
        {
            double x,y,z;
            sscanf(line+1,"%lf %lf %lf", &x, &y, &z);
            surfels.insert_vertex(tPoint3<T>(x,y,z));
        }

        else if (strncmp(line, "vn ", 3) == 0)
        {
            double x,y,z;
            sscanf(line+2,"%lf %lf %lf", &x, &y, &z);
            surfels.insert_normal(tVector3<T>(x,y,z));
        }

        else if (strncmp(line, "vc ", 3) == 0)
        {
            unsigned r,g,b;
            sscanf(line+2,"%u %u %u", &r, &g, &b);
            surfels.insert_color(ColorRgb(r/255.0, g/255.0, b/255.0));
        }
    }
}
template void read_obj(const char* name, tsurfel_set<float>& surfels);
template void read_obj(const char* name, tsurfel_set<double>& surfels);

template<class T>
void write_obj(const char* name, const tsurfel_set<T>& surfels)
{
    afree<FILE*> f(xfopen(name, "w"), fclose);

	write_obj(f, surfels);
}

template<class T>
void write_obj(FILE* f, const tsurfel_set<T>& surfels)
{
    unsigned N = surfels.size();
    for (unsigned i = 0; i < N; ++i)
    {
        const tPoint3<T>& v = surfels.vertex(i);
        fprintf(f, "v %f %f %f\n", v.x(), v.y(), v.z());

        if (surfels.has_normals())
        {
            const tVector3<T>& n = surfels.normal(i);
            fprintf(f, "vn %f %f %f\n", n.x(), n.y(), n.z());
        }

        if (surfels.has_colors())
        {
            const ColorRgb& c = surfels.vertex_color(i);
            fprintf(f, "vc %d %d %d\n",  int(c.r()*255.0), int(c.g()*255.0), int(c.b()*255.0));
        }
    }
}

template <class T>
void write_obj(const char* name, const tsurfelset_view<T>& surfels)
{
    afree<FILE*> f(xfopen(name, "w"), fclose);
    unsigned N = surfels.size();

    for (unsigned i = 0; i < N; ++i)
    {
        const tPoint3<T>& v = surfels.vertex(i);
        fprintf(f, "v %f %f %f\n", v.x(), v.y(), v.z());

        if (surfels.has_normals())
        {
            const tVector3<T>& n = surfels.normal(i);
            fprintf(f, "vn %f %f %f\n", n.x(), n.y(), n.z());
        }

        if (surfels.has_colors())
        {
            const ColorRgb& c = surfels.vertex_color(i);
            fprintf(f, "vc %d %d %d\n",  int(c.r()*255.0), int(c.g()*255.0), int(c.b()*255.0));
        }
    }
}
template void write_obj(const char* name, const tsurfelset_view<float>& surfels);
template void write_obj(const char* name, const tsurfelset_view<double>& surfels);



void write_obj(const char* name, const surfel_hierarchy_view& surfels)
{
    afree<FILE*> f(xfopen(name, "w"), fclose);
    unsigned N = surfels.size();

    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& v = surfels.vertex(i);
        fprintf(f, "v %f %f %f\n", v.x(), v.y(), v.z());

        if (surfels.has_normals())
        {
            const Vector3& n = surfels.normal(i);
            fprintf(f, "vn %f %f %f\n", n.x(), n.y(), n.z());
        }

        if (surfels.has_colors())
        {
            const ColorRgb& c = surfels.vertex_color(i);
            fprintf(f, "vc %d %d %d\n",  int(c.r()*255.0), int(c.g()*255.0), int(c.b()*255.0));
        }
    }
}

/*
 * My file format:
 *   extension .sst
 *   Format:
 *     VERSION:
 *     1  Cookie: 'SST '
 *     1  int Version - File version number
 *     1  unsigned N - number of points
 *     2  unsigned chunks - can have multiple chunks with total of N points
 *     2  unsigned K - number of points in the current chunk
 *     1  bool(byte) normals - flag that is true if have normals
 *     1  bool(byte) colors - flag that is true if have colors
 *     1  bool(byte) radius - flag that is true if have radius
 *     1  N Vertices, normals, colors and radiuses
 *
 *    Ver
 */
static const int BinPointsVersion=3;
static const unsigned BinPointsCookie='SST ';

bool read_bin_points_header(FILE* f, binpoints_header& header)
{
    fseek(f, 0, SEEK_SET);
    
    unsigned cookie;
    read_unsigned(&cookie, f);
    if (cookie != BinPointsCookie)
    {
        printf("read_bin_points_header: unknown file\n");
        return false;
    }
    read_int(&header.version, f);
    read_unsigned(&header.N, f);

    if (header.version== 1)
    {
        header.chunks=1;
    }
    else if (header.version > 1)
    {
        read_unsigned(&header.chunks, f);
    }
    read_bool(&header.hn, f);
    read_bool(&header.hc, f);
    read_bool(&header.hr, f);

    if (header.version > 2)
    {
        read_int(&header.format, f);
    }
    else
    {
        header.format = BP_FORMAT_DOUBLE;
    }

    return true;
}

bool read_bin_points_chunk_header(FILE* f, binpoints_header& header, binpoints_chunkheader& chunk_header)
{
    if (header.version == 1)
    {
        header.chunks = 1;
        chunk_header.K = header.N;
    }
    else
    {
        read_unsigned(&chunk_header.K, f);
    }
    return true;
}

bool write_bin_points_header(FILE* f, binpoints_header& header)
{
    fseek(f, 0, SEEK_SET);
    
    write_unsigned(BinPointsCookie, f);
    write_int(header.version, f);
    write_unsigned(header.N, f);
    write_unsigned(header.chunks, f);
    write_bool(header.hn, f);
    write_bool(header.hc, f);
    write_bool(header.hr, f);
    write_int(header.format, f);

    fflush(f);

    return true;
}

bool write_bin_points_chunk_header(FILE* f, const binpoints_header& header, binpoints_chunkheader& chunk_header)
{
    write_unsigned(chunk_header.K, f);
    return true;
}


//
// Read a chunch that may be of different format
//
// N - number of elements read so far
// K - number of elements to read
// hn,hn,hr - have normals, color, radius
// format - float or double
//
template <class T>
void universal_read_cunk(FILE* f, tsurfel_set<T>& ss, int N, int K, bool hn, bool hc, bool hr, int format)
{
    switch(format)
    {
    case BP_FORMAT_FLOAT:
    {
        aptr<tPoint3<float> > tpnt(new tPoint3<float>[K]);
        aptr<tVector3<float> > tvec( new tVector3<float>[K]);
        aptr<ColorRgb> col(new ColorRgb[K]);
        aptr<float> trad(new float[K]);
        fread(tpnt, sizeof(tPoint3<float>), K, f);
        if (hn) fread(tvec, sizeof(tVector3<float>), K, f);
        if (hc) fread(col, sizeof(ColorRgb), K, f);
        if (hr) fread(trad, sizeof(float), K, f);

        for (int i = 0; i < K; ++i)
        {
            ss.set_vertex(N+i, tPoint3<T>(tpnt[i].x(), tpnt[i].y(), tpnt[i].z()));
            if (hn) ss.set_normal(N+i, tVector3<T>(tvec[i].x(), tvec[i].y(), tvec[i].z()));
            if (hc) ss.set_color(N+i, col[i]);
            if (hr) ss.set_radius(N+i, trad[i]);
        }
        break;
    }
    case BP_FORMAT_DOUBLE:
    {
        aptr<tPoint3<double> > tpnt(new tPoint3<double>[K]);
        aptr<tVector3<double> > tvec( new tVector3<double>[K]);
        aptr<ColorRgb> col(new ColorRgb[K]);
        aptr<double> trad(new double[K]);
        fread(tpnt, sizeof(tPoint3<double>), K, f);
        if (hn) fread(tvec, sizeof(tVector3<double>), K, f);
        if (hc) fread(col, sizeof(ColorRgb), K, f);
        if (hr) fread(trad, sizeof(double), K, f);

        for (int i = 0; i < K; ++i)
        {
            ss.set_vertex(N+i, tPoint3<T>(tpnt[i].x(), tpnt[i].y(), tpnt[i].z()));
            if (hn) ss.set_normal(N+i, tVector3<T>(tvec[i].x(), tvec[i].y(), tvec[i].z()));
            if (hc) ss.set_color(N+i, col[i]);
            if (hr) ss.set_radius(N+i, trad[i]);
        }
        break;
    } // case
    } //switch
}

template void universal_read_cunk(FILE* f, tsurfel_set<float>& ss, int N, int K, bool hn, bool hc, bool hr, int format);
template void universal_read_cunk(FILE* f, tsurfel_set<double>& ss, int N, int K, bool hn, bool hc, bool hr, int format);

template <class T>
void read_bin_points(FILE* f, tsurfel_set<T>& surfels)
{
    bool hn,hc,hr;
    unsigned N,chunks;
    int currentN = 0; // number of points read so far
    binpoints_header header;
    if (!read_bin_points_header(f, header)) return;

    N=header.N;
    chunks = header.chunks;
    hn = header.hn;
    hc = header.hc;
    hr = header.hr;

    surfels.resize(N, hn ? N : 0, hc ? N : 0, hr ? N : 0);
    int my_format = sizeof(T) == 4 ? BP_FORMAT_FLOAT : BP_FORMAT_DOUBLE; // Q&D HACK

    while (chunks--)
    {
        binpoints_chunkheader chunk_header;
        read_bin_points_chunk_header(f, header, chunk_header);
        unsigned K = chunk_header.K;
        if (my_format == header.format)
        {
            fread(&(surfels.vertex(currentN)[0]), sizeof(typename tsurfel_set<T>::vertex_list::value_type), K, f);
            if (hn) fread(&(surfels.normal(currentN)[0]), sizeof(typename tsurfel_set<T>::normal_list::value_type), K, f);
            if (hc) fread(&(surfels.vertex_color(currentN)[0]), sizeof(typename tsurfel_set<T>::color_list::value_type), K, f);
            if (hr) fread(&(surfels.radiuses()[currentN]), sizeof(typename tsurfel_set<T>::radius_list::value_type), K, f);
        }
        else
        {
            universal_read_cunk(f, surfels, currentN, K, hn, hc, hr, header.format);
        }
        currentN += K;
    }
}

template <class T>
void read_bin_points(const char* name, tsurfel_set<T>& surfels)
{
    afree<FILE*> f(fopen(name, "rb"), fclose);
    if (f == 0)
    {
        printf("read_bin_points: Failed to open file for read\n");
        return;
    }

	read_bin_points(f, surfels);
}

template void read_bin_points(const char* name, tsurfel_set<float>& surfels);
template void read_bin_points(const char* name, tsurfel_set<double>& surfels);

template<class T>
void write_bin_points(FILE* f, const tsurfel_set<T>& surfels)
{
    unsigned N = surfels.size();
    bool hn = surfels.has_normals();
    bool hc = surfels.has_colors();
    bool hr = surfels.has_radius();
    binpoints_header header;
    header.version = BinPointsVersion;
    header.N = N;
    header.chunks = 1;
    header.hn = hn;
    header.hc = hc;
    header.hr = hr;
    header.format = sizeof(T)==4 ? BP_FORMAT_FLOAT : BP_FORMAT_DOUBLE; // HACK so I don't have to specialize for float and double

    write_bin_points_header(f, header);

    binpoints_chunkheader chunk_header;
    chunk_header.K = N;
    write_bin_points_chunk_header(f, header, chunk_header);

    fwrite(&(surfels.vertices()[0]), sizeof(typename tsurfel_set<T>::vertex_list::value_type), N, f);
    if (hn) fwrite(&(surfels.normals()[0]), sizeof(typename tsurfel_set<T>::normal_list::value_type), N, f);
    if (hc) fwrite(&(surfels.vertex_colors()[0]), sizeof(typename tsurfel_set<T>::color_list::value_type), N, f);
    if (hr) fwrite(&(surfels.radiuses()[0]), sizeof(typename tsurfel_set<T>::radius_list::value_type), N, f);
}

template void write_bin_points(FILE* f, const tsurfel_set<float>& surfels);
template void write_bin_points(FILE* f, const tsurfel_set<double>& surfels);

template <class T>
void append_bin_points(FILE* f, const tsurfel_set<T>& surfels)
{
    unsigned N = surfels.size();
    bool hn = surfels.has_normals();
    bool hc = surfels.has_colors();
    bool hr = surfels.has_radius();

// bugbug check format
    binpoints_header header;
    if (!read_bin_points_header(f, header)) return;
    if (header.version == 1)
    {
        printf("Cannot append to old files\n");
        return;
    }
    if ((hn != header.hn) || (hc != header.hc) || (hr != header.hr))
    {
        printf("append_bin_points: incompatible file, i.e. missing normals or colors\n");
        return;
    }
    header.N += N;
    ++header.chunks;
    write_bin_points_header(f, header);

    fseek(f, 0, SEEK_END);

    binpoints_chunkheader chunk_header;
    chunk_header.K = N;
    write_bin_points_chunk_header(f, header, chunk_header);

    fwrite(&(surfels.vertices()[0]), sizeof(typename tsurfel_set<T>::vertex_list::value_type), N, f);
    if (hn) fwrite(&(surfels.normals()[0]), sizeof(typename tsurfel_set<T>::normal_list::value_type), N, f);
    if (hc) fwrite(&(surfels.vertex_colors()[0]), sizeof(typename tsurfel_set<T>::color_list::value_type), N, f);
    if (hr) fwrite(&(surfels.radiuses()[0]), sizeof(typename tsurfel_set<T>::radius_list::value_type), N, f);
}

template void append_bin_points(FILE* f, const tsurfel_set<float>& surfels);
template void append_bin_points(FILE* f, const tsurfel_set<double>& surfels);

template <class T>
void append_bin_points(const char* name, const tsurfel_set<T>& surfels)
{
    afree<FILE*> f(fopen(name, "rb+"), fclose);
    if (f == 0)
    {
        printf("read_bin_points: Failed to open file for read\n");
        return;
    }

	append_bin_points(f, surfels);
}

template void append_bin_points(const char* name, const tsurfel_set<float>& surfels);
template void append_bin_points(const char* name, const tsurfel_set<double>& surfels);


template <class T>
void write_bin_points(const char* name, const tsurfel_set<T>& surfels)
{
    afree<FILE*> f(fopen(name, "wb"), fclose);
    if (f == 0)
    {
        printf("write_bin_points: Failed to open file for write\n");
        return;
    }

	write_bin_points(f, surfels);
}

template void write_bin_points(const char* name, const tsurfel_set<float>& surfels);
template void write_bin_points(const char* name, const tsurfel_set<double>& surfels);


void print(const char* prefix, const surfel_set& surfels)
{
    unsigned n_vertices = surfels.size();

    for (unsigned i = 0; i < n_vertices; ++i)
    {
        const Point3& p = surfels.vertex(i);
        printf("%s [ %1.20g %1.20g %1.20g ]\n", prefix, p[0], p[1], p[2]);
    }
}

void print(const char* prefix, const surfelset_view& sv)
{
    unsigned  N = sv.size();
    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& p = sv.vertex(i);
        printf("%s [ %1.20g %1.20g %1.20g ]\n", prefix, p[0], p[1], p[2]);
    }
}

void print(const char* prefix, const surfel_hierarchy& hierarchy)
{
    unsigned num_levels = hierarchy.num_levels();
    for (unsigned level = 0; level < num_levels; ++level)
    {
        print(prefix, hierarchy.level(level));
    }
}

void print(const char* prefix, const surfel_hierarchy_view& hv)
{
    unsigned  N = hv.size();
    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& p = hv.vertex(i);
        printf("%s [ %1.20g %1.20g %1.20g ]\n", prefix, p[0], p[1], p[2]);
    }
}

void print(const char* prefix, const point3_set& p3s)
{
    unsigned  N = p3s.size();
    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& p = p3s[i];
        printf("%s [ %1.20g %1.20g %1.20g ]\n", prefix, p[0], p[1], p[2]);
    }
}

GTB_END_NAMESPACE
