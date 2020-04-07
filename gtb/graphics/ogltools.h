
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


#ifndef __OGLTOOLS_H
#define __OGLTOOLS_H

#include <gtb/common.hpp>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <gtb/graphics/point3.hpp>

GTB_BEGIN_NAMESPACE

extern unsigned char ctable[][3];
extern const char *cnames[];
extern int ncolors;

void print_gl_errors();

template <class T>
struct OpenGL
{
    static void Vertex(T, T);
    static void Vertex2v(const T*);
    static void Vertex(T, T, T);
    static void Vertex3v(const T*);
    static void Vertex(T, T, T, T);
    static void Vertex4v(const T*);

    static void Color(T, T, T);
    static void Color3v(const T*);
    static void Color(T, T, T, T);
    static void Color4v(const T*);

    static void Normal(T, T, T);
    static void Normal3v(const T*);
private:
    OpenGL() {};
    OpenGL(const OpenGL&) {};
};

struct CInfo
{
    typedef unsigned char value_type;
    value_type r,g,b;

    CInfo() {}
    CInfo(const CInfo& rhs) : r(rhs.r), g(rhs.g), b(rhs.b) {}
    CInfo(value_type R, value_type G, value_type B) : r(R), g(G), b(B) {}
    CInfo(double R, double G, double B) :
	    r((unsigned char) (R*255.0)),
	    g((unsigned char) (G*255.0)),
	    b((unsigned char) (B*255.0)) {}
    CInfo(value_type* c) : r(c[0]), g(c[1]), b(c[2]) {}

    CInfo Mix( const CInfo& rhs, double ratio) const
    {
        double ratio1 = 1-ratio;
        return CInfo(
            (value_type)(r * ratio + rhs.r * ratio1),
            (value_type)(g * ratio + rhs.g * ratio1),
            (value_type)(b * ratio + rhs.b * ratio1));
    }

    operator value_type*() { return &r; }
    void load() const { glColor3ub(r,g,b); }
    void ogl() const { glColor3ub(r,g,b); }
    void ogl(value_type a) const { glColor4ub(r,g,b,a); }

    void load_as_clearcolor() const
    {
        glClearColor(r/255.0, g/255.0, b/255.0, 0);
    }
};



void WriteOGLBuffer(const char* name);
void WriteOGLZBuffer(const char* name);
int GetColorIndex(const char* name);
void WriteFB();
Point3 ProjectCurrent(Point3& p);
void CheckGLE();
void print_gl_info();
int next_color(int color);
int prev_color(int color);

//
// read into memory / write back ogl buffers
//
void* get_depth_buffer(GLint dims[4]);
void* get_color_buffer(GLint dims[4]);
void put_depth_buffer(void* frame, GLint dims[4]);
void put_color_buffer(void* frame, GLint dims[4]);

inline CInfo GetCInfo(int idx) { return CInfo(ctable[idx]); }

inline CInfo MixColors(int idx1, int idx2, double ratio)
{
    CInfo mixed = GetCInfo(idx1).Mix(GetCInfo(idx2), ratio);
    return mixed;
}

/*
 * 0..0.5 blue = 1 , red varies linearly
 * 0.5..1 reversed
 */
inline CInfo MixRedBlue(double ratio)
{
    CInfo c;
    c.r = ratio < 0.5 ? (CInfo::value_type)(255*2*ratio) : 255;
    c.g = 0;
    c.b = ratio < 0.5 ? 255 : (CInfo::value_type)(255*(2-ratio*2));
    return c;
}

/*
 * Mix colors to create some rainbow (colormaps)
 * x - in the range [0,1]
 * N - number of colors in the color table to mix
 * colors - indices of colors to mix
 * ranges - of colors to mix
 */
CInfo MixRainbow(double x, int N, int* colors, double* ranges);
CInfo MixRainbow(double x, int rainbowindex);
void RainbowMap(CInfo* image, int w, int h, int rainbowindex);

// The rainbows info
extern int *rainbows[]; // The colors
extern double *rbrs[];  // A uniform distribution for the above
extern int rainbowsizes[]; // Their number of colors in each

extern int rainbow1[]; // smooth blue to white:      b c g y r m w
extern int rainbow2[]; // smooth white to red        w m b c g y r
extern int rainbow3[]; // non smooth white to red    w b c g m y r
extern int rainbow4[]; // half colors
extern double rainbowranges[];



/*
 * Rainbow mix
 * Input a value in the range[0,1]
 * Output a color
 */
inline CInfo MixRainbow0(double q)
{
#define V1  1e-6
#define V2  0.1
#define V3  0.33
#define V4  0.66
#define V5  1.0

    double r=0,g=0,b=0;

    if (q < V1) {r = 1; g=1; b=1; }
    else if (q < V2) r = 1;
    else if (q < V3) r = 1-(q-V3)/(V3-V2);

    if ((q>V1) && (q<V2)) g = (q-V1)/(V2-V1);
    if ((q>V2) && (q<V4)) g = 1;
    else if (q>=V4) g = 1.0 - (q-V4)/(V5-V4); 

    if (q>V3) b = (q-V3)/(V4-V3);
    else if (q>=V4) b = 1;

    return CInfo(r,g,b);
#undef V1
#undef V2
#undef V3
#undef V4
#undef V5
}

inline void LoadMixColors(int idx1, int idx2, double ratio)
{
    MixColors(idx1, idx2, ratio).load();
}

inline void LoadGLColor(int index, unsigned char alpha = 255)
{
    if (alpha != 255)
    {
        unsigned char* c = ctable[index % ncolors];
        glColor4ub(c[0], c[1], c[2], alpha);
    }
    else
    {
        glColor3ubv(ctable[index % ncolors]);
    }
}

inline void LoadAsBackgroundColor(int index)
{
    GetCInfo(index).load_as_clearcolor();
}

void bitmap_output(double x, double y, char *string, void *font);

enum COLORNAMES {
    Snow=0,
    GhostWhite,
    WhiteSmoke,
    Gainsboro,
    FloralWhite,
    OldLace,
    Linen,
    AntiqueWhite,
    PapayaWhip,
    BlanchedAlmond,
    Bisque,
    PeachPuff,
    NavajoWhite,
    Moccasin,
    Cornsilk,
    Ivory,
    LemonChiffon,
    Seashell,
    Honeydew,
    MintCream,
    Azure,
    AliceBlue,
    lavender,
    LavenderBlush,
    MistyRose,
    White,
    Black,
    DarkSlateGray,
    DimGrey,
    SlateGrey,
    LightSlateGray,
    Grey,
    LightGray,
    MidnightBlue,
    NavyBlue,
    CornflowerBlue,
    DarkSlateBlue,
    SlateBlue,
    MediumSlateBlue,
    LightSlateBlue,
    MediumBlue,
    RoyalBlue,
    Blue,
    DodgerBlue,
    DeepSkyBlue,
    SkyBlue,
    LightSkyBlue,
    SteelBlue,
    LightSteelBlue,
    LightBlue,
    PowderBlue,
    PaleTurquoise,
    DarkTurquoise,
    MediumTurquoise,
    Turquoise,
    Cyan,
    LightCyan,
    CadetBlue,
    MediumAquamarine,
    Aquamarine,
    DarkGreen,
    DarkOliveGreen,
    DarkSeaGreen,
    SeaGreen,
    MediumSeaGreen,
    LightSeaGreen,
    PaleGreen,
    SpringGreen,
    LawnGreen,
    Green,
    MonsterGreen,
    Chartreuse,
    MedSpringGreen,
    GreenYellow,
    LimeGreen,
    YellowGreen,
    ForestGreen,
    OliveDrab,
    DarkKhaki,
    PaleGoldenrod,
    LtGoldenrodYello,
    LightYellow,
    Yellow,
    Gold,
    LightGoldenrod,
    Goldenrod,
    DarkGoldenrod,
    RosyBrown,
    IndianRed,
    SaddleBrown,
    Sienna,
    Peru,
    Burlywood,
    Beige,
    Wheat,
    SandyBrown,
    Tan,
    Chocolate,
    Firebrick,
    Brown,
    DarkSalmon,
    Salmon,
    LightSalmon,
    Orange,
    DarkOrange,
    Coral,
    LightCoral,
    Tomato,
    OrangeRed,
    Red,
    HotPink,
    DeepPink,
    Pink,
    LightPink,
    PaleVioletRed,
    Maroon,
    MediumVioletRed,
    VioletRed,
    Magenta,
    Violet,
    Plum,
    Orchid,
    MediumOrchid,
    DarkOrchid,
    DarkViolet,
    BlueViolet,
    Purple,
    MediumPurple,
    Thistle,
    Snow1,
    Snow2,
    Snow3,
    Snow4,
    Seashell1,
    Seashell2,
    Seashell3,
    Seashell4,
    AntiqueWhite1,
    AntiqueWhite2,
    AntiqueWhite3,
    AntiqueWhite4,
    Bisque1,
    Bisque2,
    Bisque3,
    Bisque4,
    PeachPuff1,
    PeachPuff2,
    PeachPuff3,
    PeachPuff4,
    NavajoWhite1,
    NavajoWhite2,
    NavajoWhite3,
    NavajoWhite4,
    LemonChiffon1,
    LemonChiffon2,
    LemonChiffon3,
    LemonChiffon4,
    Cornsilk1,
    Cornsilk2,
    Cornsilk3,
    Cornsilk4,
    Ivory1,
    Ivory2,
    Ivory3,
    Ivory4,
    Honeydew1,
    Honeydew2,
    Honeydew3,
    Honeydew4,
    LavenderBlush1,
    LavenderBlush2,
    LavenderBlush3,
    LavenderBlush4,
    MistyRose1,
    MistyRose2,
    MistyRose3,
    MistyRose4,
    Azure1,
    Azure2,
    Azure3,
    Azure4,
    SlateBlue1,
    SlateBlue2,
    SlateBlue3,
    SlateBlue4,
    RoyalBlue1,
    RoyalBlue2,
    RoyalBlue3,
    RoyalBlue4,
    Blue1,
    Blue2,
    Blue3,
    Blue4,
    DodgerBlue1,
    DodgerBlue2,
    DodgerBlue3,
    DodgerBlue4,
    SteelBlue1,
    SteelBlue2,
    SteelBlue3,
    SteelBlue4,
    DeepSkyBlue1,
    DeepSkyBlue2,
    DeepSkyBlue3,
    DeepSkyBlue4,
    SkyBlue1,
    SkyBlue2,
    SkyBlue3,
    SkyBlue4,
    LightSkyBlue1,
    LightSkyBlue2,
    LightSkyBlue3,
    LightSkyBlue4,
    SlateGray1,
    SlateGray2,
    SlateGray3,
    SlateGray4,
    LightSteelBlue1,
    LightSteelBlue2,
    LightSteelBlue3,
    LightSteelBlue4,
    LightBlue1,
    LightBlue2,
    LightBlue3,
    LightBlue4,
    LightCyan1,
    LightCyan2,
    LightCyan3,
    LightCyan4,
    PaleTurquoise1,
    PaleTurquoise2,
    PaleTurquoise3,
    PaleTurquoise4,
    CadetBlue1,
    CadetBlue2,
    CadetBlue3,
    CadetBlue4,
    Turquoise1,
    Turquoise2,
    Turquoise3,
    Turquoise4,
    Cyan1,
    Cyan2,
    Cyan3,
    Cyan4,
    DarkSlateGray1,
    DarkSlateGray2,
    DarkSlateGray3,
    DarkSlateGray4,
    Aquamarine1,
    Aquamarine2,
    Aquamarine3,
    Aquamarine4,
    DarkSeaGreen1,
    DarkSeaGreen2,
    DarkSeaGreen3,
    DarkSeaGreen4,
    SeaGreen1,
    SeaGreen2,
    SeaGreen3,
    SeaGreen4,
    PaleGreen1,
    PaleGreen2,
    PaleGreen3,
    PaleGreen4,
    SpringGreen1,
    SpringGreen2,
    SpringGreen3,
    SpringGreen4,
    Green1,
    Green2,
    Green3,
    Green4,
    Chartreuse1,
    Chartreuse2,
    Chartreuse3,
    Chartreuse4,
    OliveDrab1,
    OliveDrab2,
    OliveDrab3,
    OliveDrab4,
    DarkOliveGreen1,
    DarkOliveGreen2,
    DarkOliveGreen3,
    DarkOliveGreen4,
    Khaki1,
    Khaki2,
    Khaki3,
    Khaki4,
    LightGoldenrod1,
    LightGoldenrod2,
    LightGoldenrod3,
    LightGoldenrod4,
    LightYellow1,
    LightYellow2,
    LightYellow3,
    LightYellow4,
    Yellow1,
    Yellow2,
    Yellow3,
    Yellow4,
    Gold1,
    Gold2,
    Gold3,
    Gold4,
    Goldenrod1,
    Goldenrod2,
    Goldenrod3,
    Goldenrod4,
    DarkGoldenrod1,
    DarkGoldenrod2,
    DarkGoldenrod3,
    DarkGoldenrod4,
    RosyBrown1,
    RosyBrown2,
    RosyBrown3,
    RosyBrown4,
    IndianRed1,
    IndianRed2,
    IndianRed3,
    IndianRed4,
    Sienna1,
    Sienna2,
    Sienna3,
    Sienna4,
    Burlywood1,
    Burlywood2,
    Burlywood3,
    Burlywood4,
    Wheat1,
    Wheat2,
    Wheat3,
    Wheat4,
    Tan1,
    Tan2,
    Tan3,
    Tan4,
    Chocolate1,
    Chocolate2,
    Chocolate3,
    Chocolate4,
    Firebrick1,
    Firebrick2,
    Firebrick3,
    Firebrick4,
    Brown1,
    Brown2,
    Brown3,
    Brown4,
    Salmon1,
    Salmon2,
    Salmon3,
    Salmon4,
    LightSalmon1,
    LightSalmon2,
    LightSalmon3,
    LightSalmon4,
    Orange1,
    Orange2,
    Orange3,
    Orange4,
    DarkOrange1,
    DarkOrange2,
    DarkOrange3,
    DarkOrange4,
    Coral1,
    Coral2,
    Coral3,
    Coral4,
    Tomato1,
    Tomato2,
    Tomato3,
    Tomato4,
    OrangeRed1,
    OrangeRed2,
    OrangeRed3,
    OrangeRed4,
    Red1,
    Red2,
    Red3,
    Red4,
    DeepPink1,
    DeepPink2,
    DeepPink3,
    DeepPink4,
    HotPink1,
    HotPink2,
    HotPink3,
    HotPink4,
    Pink1,
    Pink2,
    Pink3,
    Pink4,
    LightPink1,
    LightPink2,
    LightPink3,
    LightPink4,
    PaleVioletRed1,
    PaleVioletRed2,
    PaleVioletRed3,
    PaleVioletRed4,
    Maroon1,
    Maroon2,
    Maroon3,
    Maroon4,
    VioletRed1,
    VioletRed2,
    VioletRed3,
    VioletRed4,
    Magenta1,
    Magenta2,
    Magenta3,
    Magenta4,
    Orchid1,
    Orchid2,
    Orchid3,
    Orchid4,
    Plum1,
    Plum2,
    Plum3,
    Plum4,
    MediumOrchid1,
    MediumOrchid2,
    MediumOrchid3,
    MediumOrchid4,
    DarkOrchid1,
    DarkOrchid2,
    DarkOrchid3,
    DarkOrchid4,
    Purple1,
    Purple2,
    Purple3,
    Purple4,
    MediumPurple1,
    MediumPurple2,
    MediumPurple3,
    MediumPurple4,
    Thistle1,
    Thistle2,
    Thistle3,
    Thistle4,
    Grey11,
    Grey21,
    Grey31,
    Grey41,
    Grey51,
    Grey61,
    Grey71,
    Grey81,
    Grey91,
    DarkGrey,
    DarkBlue,
    DarkCyan,
    DarkMagenta,
    DarkRed,
    LightGreen
};

extern int ncolors; // Number of colors in the list

struct MatrialInfo
{
    const char* name;
    float amb[4];
    float diff[4];
    float spec[4];
    float shininess;
};

extern MatrialInfo MaterialsList[];
extern const int nMaterials;
void load_material(int matidx, GLenum face);
void reset_material(GLenum face);

enum MATERIALNAMES {
    mtBrass,
    mtBronze,
    mtPolishedBronze,
    mtChrome,
    mtCopper,
    mtPolishedCopper,
    mtGold,
    mtPolishedGold,
    mtPewter,
    mtSilver,
    mtPolishedSilver,
    mtEmerald,
    mtJade,
    mtObsidian,
    mtPearl,
    mtRuby,
    mtTurquoise,
    mtBlackPlastic,
    mtBlackRubber,
    mt2emerald,
    mt2jade,
    mt2obsidian,
    mt2pearl,
    mt2ruby,
    mt2turquoise,
    mt2brass,
    mt2bronze,
    mt2chrome,
    mt2copper,
    mt2gold,
    mt2silver,
    mt2blackplastic,
    mt2cyanplastic,
    mt2greenplastic,
    mt2redplastic,
    mt2whiteplastic,
    mt2yellowplastic,
    mt2blackrubber,
    mt2cyanrubber,
    mt2greenrubber,
    mt2redrubber,
    mt2whiterubber,
    mt2yellowrubber,
};


GTB_END_NAMESPACE

#include <gtb/graphics/ogltools.inl>

#endif // __OGLTOOLS_H
