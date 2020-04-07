
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
#include <gtb/graphics/ogltools.h>
#include <gtb/graphics/camera.hpp>
#include <gtb/graphics/line3.hpp>
#include <gtb/error/error.hpp>
#include <gtb/io/io.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/camera.ipp>
#undef inline
#endif


using namespace std;

GTB_BEGIN_NAMESPACE


static const double DEFAULT_X_FOV = M_PI / 6.0;
static const double DEFAULT_Y_FOV = M_PI / 6.0;

enum {
	NEAR_PLANE,
	FAR_PLANE,
	LEFT_PLANE,
	RIGHT_PLANE,
	TOP_PLANE,
	BOTTOM_PLANE
};


template <class T>
tCamera<T>::tCamera()
{
	_row = 0;
	_col = 0;
	_n_rows = 1;
	_n_cols = 1;
	_type = PERSPECTIVE;

	reset(Point3(0.0, 0.0, 0.0),
	      Vector3(0.0, 0.0, -1.0),
	      Vector3(0.0, 1.0, 0.0),
	      (T)DEFAULT_X_FOV,
	      (T)DEFAULT_Y_FOV,
	      (T)0.001,
	      (T)1000.0);
}


template <class T>
tCamera<T>::tCamera(const Point3 &c_origin,
	       const Vector3 &c_towards,
	       const Vector3 &c_up,
	       value_type c_x_fov,
	       value_type c_y_fov,
	       value_type c_near_distance,
	       value_type c_far_distance)
{
	_row = 0;
	_col = 0;
	_n_rows = 1;
	_n_cols = 1;
	_type = PERSPECTIVE;

	reset(c_origin,
	      c_towards,
	      c_up,
	      c_x_fov,
	      c_y_fov,
	      c_near_distance,
	      c_far_distance);
}


template <class T>
tCamera<T>::tCamera(const tCamera &cam)
{
	*this = cam;
}


template <class T>
tCamera<T>::tCamera(const tCamera &cam,
	       unsigned i,
	       unsigned j,
	       unsigned n_rows,
	       unsigned n_cols)
{
	_row = i;
	_col = j;
	_n_rows = n_rows;
	_n_cols = n_cols;
	_type = cam._type;

	assert(_row < _n_rows);
	assert(_col < _n_cols);
	assert((_n_rows > 0) && ((_n_rows == 1) || ((_n_rows % 2) == 0)));
	assert((_n_cols > 0) && ((_n_cols == 1) || ((_n_cols % 2) == 0)));

	reset(cam.origin(),
	      cam.towards(),
	      cam.up(),
	      cam.x_fov(),
	      cam.y_fov(),
	      cam.near_distance(),
	      cam.far_distance());
}


template <class T>
tCamera<T> &tCamera<T>::operator=(const tCamera &cam)
{
	if (&cam != this) {
		_cs = cam._cs;
		_x_fov = cam._x_fov;
		_y_fov = cam._y_fov;
		_near_distance = cam._near_distance;
		_far_distance = cam._far_distance;
		for (unsigned i = 0; i < 6; i++) {
			_planes[i] = cam._planes[i];
		}
		for (unsigned i = 0; i < 8; i++) {
			_corners[i] = cam._corners[i];
		}
		_row = cam._row;
		_col = cam._col;
		_n_rows = cam._n_rows;
		_n_cols = cam._n_cols;
		_type = cam._type;
	}
	return *this;
}


template <class T>
void tCamera<T>::reset(const Point3 &arg_origin,
		   const Vector3 &arg_towards,
		   const Vector3 &arg_up)
{
	reset(arg_origin,
	      arg_towards,
	      arg_up,
	      x_fov(),
	      y_fov(),
	      near_distance(),
	      far_distance());
}


template <class T>
void tCamera<T>::reset(const Point3 &arg_origin,
		   const Vector3 &arg_towards,
		   const Vector3 &arg_up,
		   value_type arg_x_fov,
		   value_type arg_y_fov,
		   value_type arg_near_distance,
		   value_type arg_far_distance)
{
	Vector3 z = -arg_towards;
	z.normalize();
	Vector3 x = arg_up.cross(z);
	x.normalize();
	Vector3 y = z.cross(x);
	_cs.reset(arg_origin, x, y, z);

	_x_fov = arg_x_fov;
	_y_fov = arg_y_fov;
	_near_distance = arg_near_distance;
	_far_distance = arg_far_distance;
	
	assert(_row < _n_rows);
	assert(_col < _n_cols);
	assert((_n_rows > 0) && ((_n_rows == 1) || ((_n_rows % 2) == 0)));
	assert((_n_cols > 0) && ((_n_cols == 1) || ((_n_cols % 2) == 0)));
	Point3 O = origin();
	Point3 O1 = O + (near_distance() * towards());
	Point3 O2 = O + (far_distance() * towards());
	Point3 *A = &_corners[0];
	Point3 *B = &_corners[1];
	Point3 *C = &_corners[2];
	Point3 *D = &_corners[3];
	Point3 *E = &_corners[4];
	Point3 *F = &_corners[5];
	Point3 *G = &_corners[6];
	Point3 *H = &_corners[7];

	value_type x1 = near_distance() * tan(x_fov());
	value_type y1 = near_distance() * tan(y_fov());

	value_type x2 = far_distance() * tan(x_fov());
	value_type y2 = far_distance() * tan(y_fov());

	Vector3 r1 = right() * (T)(2 * x1 / (T)_n_cols);
	Vector3 d1 = down() * (T)(2 * y1 / (T)_n_rows);

	Vector3 r2 = right() * (T)(2 * x2 / (T)_n_cols);
	Vector3 d2 = down() * (T)(2 * y2 / (T)_n_rows);

	*B = O1 + (x1 * left()) + (y1 * up()) + ((value_type)_col * r1) + ((value_type)_row * d1);
	*C = *B + r1;
	*A = *B + d1;
	*D = *B + r1 + d1;

	*F = O2 + (x2 * left()) + (y2 * up()) + ((value_type)_col * r2) + ((value_type)_row * d2);
	*G = *F + r2;
	*E = *F + d2;
	*H = *F + r2 + d2;

	_planes[NEAR_PLANE] = Plane(*A, *B, *C);
	_planes[FAR_PLANE] = Plane(*F, *E, *H);
	_planes[LEFT_PLANE] = Plane(*A, *E, *F);
	_planes[RIGHT_PLANE] = Plane(*C, *G, *H);
	_planes[TOP_PLANE] = Plane(*B, *F, *G);
	_planes[BOTTOM_PLANE] = Plane(*E, *A, *D);

	assert(!sees(origin()));
}


template <class T>
void tCamera<T>::setup_viewport(const Viewport &viewport) const
{
//  	int x = (int) ((value_type) _col / (value_type) _n_cols *
//  		       (value_type) viewport.width());
//  	int y = (int) ((value_type) (_n_rows - 1 - _row) / (value_type) _n_rows *
//  		       (value_type) viewport.height());
//  	int w = viewport.width() / _n_cols;
//  	int h = viewport.height() / _n_rows;
//  	glViewport(x, y, w, h);
  	glViewport(viewport.x_min(), viewport.y_min(),
		   viewport.width(), viewport.height());
}


template <class T>
void tCamera<T>::setup_frustum() const
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	GLdouble ymax = _near_distance * tan(_y_fov);
	GLdouble ymin = -ymax;
	GLdouble xmax = _near_distance * tan(_x_fov);
	GLdouble xmin = -xmax;

	GLdouble w = xmax - xmin;
	GLdouble h = ymax - ymin;
	GLdouble dw = w / _n_cols;
	GLdouble dh = h / _n_rows;

	GLdouble l = xmin + _col * dw;
	GLdouble r = l + dw;
	GLdouble t = ymax - _row * dh;
	GLdouble b = t - dh;

	glFrustum(l, r, b, t, _near_distance, _far_distance);
}


template <class T>
void tCamera<T>::load_headlights(const vector<Light> &headlights) const
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	for (unsigned i = 0; i < headlights.size(); i++) {
		headlights[i].load();
	}

}


template <class T>
void tCamera<T>::setup_look_at() const
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	Point3 look_at = origin() + near_distance() * towards();
	gluLookAt(origin().x(), origin().y(), origin().z(),
		  look_at.x(), look_at.y(), look_at.z(),
		  up().x(), up().y(), up().z());
}


template <class T>
void tCamera<T>::load(const Viewport &viewport) const
{
	switch (_type) {
	case PERSPECTIVE:
		setup_viewport(viewport);
		setup_frustum();
		setup_look_at();
		break;
	case ORTHOGRAPHIC:
		viewport.load();

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(viewport.x_min(), viewport.width(),
			   viewport.y_min(), viewport.height());

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		break;
	default:
		GTB_ERROR("invalid type");
		break;
	}
}

template <class T>
void tCamera<T>::gl_load(const Viewport &viewport) const
{
	switch (_type) {
	case PERSPECTIVE:
		setup_viewport(viewport);
		setup_frustum();
		setup_look_at();
		break;
	case ORTHOGRAPHIC:
		viewport.load();

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(viewport.x_min(), viewport.width(),
			   viewport.y_min(), viewport.height());

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		break;
	default:
		GTB_ERROR("invalid type");
		break;
	}
}


template <class T>
void tCamera<T>::load(const Viewport &viewport,
		  const vector<Light> &headlights) const
{
	// @@@ FIXME: handle orthographic case
	assert(_type == PERSPECTIVE);
	setup_viewport(viewport);
	setup_frustum();
	load_headlights(headlights);
	setup_look_at();
}


//  void tCamera<T>::render() const
//  {
//  	const Point3 &O = origin();
//  	const Point3 &A = _corners[0];
//  	const Point3 &B = _corners[1];
//  	const Point3 &C = _corners[2];
//  	const Point3 &D = _corners[3];
//  	const Point3 &E = _corners[4];
//  	const Point3 &F = _corners[5];
//  	const Point3 &G = _corners[6];
//  	const Point3 &H = _corners[7];

//  	// near
//  	glBegin(GL_LINE_LOOP);
//  	A.load();
//  	B.load();
//  	C.load();
//  	D.load();
//  	glEnd();
	
//  	// far
//  	glBegin(GL_LINE_LOOP);
//  	E.load();
//  	F.load();
//  	G.load();
//  	H.load();
//  	glEnd();

//  	// left
//  	glBegin(GL_LINE_LOOP);
//  	O.load();
//  	E.load();
//  	F.load();
//  	glEnd();

//  	// right
//  	glBegin(GL_LINE_LOOP);
//  	O.load();
//  	G.load();
//  	H.load();
//  	glEnd();

//  	// normal scale
//  	Point3 p1 = centroid(A, B, E, F);
//  	Point3 p2 = centroid(C, D, G, H);
//  	value_type l = Point3::distance(p1, p2) / 3.0;

//  	// left normal
//  	Point3 q = p1 + (l * _planes[LEFT_PLANE].normal());
//  	glBegin(GL_LINES);
//  	p1.load();
//  	q.load();
//  	glEnd();

//  //  	// right normal
//  //  	q = p2 + (l * _planes[RIGHT_PLANE].normal());
//  //  	glBegin(GL_LINES);
//  //  	p2.load();
//  //  	q.load();
//  //  	glEnd();

//  //  	// top normal
//  //  	p1 = centroid(B, C, F, G);
//  //  	q = p1 + (l * _planes[TOP_PLANE].normal());
//  //  	glBegin(GL_LINES);
//  //  	p1.load();
//  //  	q.load();
//  //  	glEnd();

//  //  	// bottom normal
//  //  	p1 = centroid(A, D, E, H);
//  //  	q = p1 + (l * _planes[BOTTOM_PLANE].normal());
//  //  	glBegin(GL_LINES);
//  //  	p1.load();
//  //  	q.load();
//  //  	glEnd();

//  //  	// near normal
//  //  	p1 = centroid(A, B, C, D);
//  //  	q = p1 + (l * _planes[NEAR_PLANE].normal());
//  //  	glBegin(GL_LINES);
//  //  	p1.load();
//  //  	q.load();
//  //  	glEnd();

//  	// far normal
//  	p1 = centroid(E, F, G, H);
//  	q = p1 + (10 * l * _planes[FAR_PLANE].normal());
//  	glBegin(GL_LINES);
//  	p1.load();
//  	q.load();
//  	glEnd();
//  }


template <class T>
void tCamera<T>::render(value_type scale) const
{
#if 1
	Point3 O = origin();
	Point3 A = O + ((_corners[0] - O) * scale);
	Point3 B = O + ((_corners[1] - O) * scale);
	Point3 C = O + ((_corners[2] - O) * scale);
	Point3 D = O + ((_corners[3] - O) * scale);

	glBegin(GL_LINE_LOOP);
	A.load();
	B.load();
	C.load();
	D.load();
	A.load();
	O.load();
	D.load();
	C.load();
	O.load();
	B.load();
	glEnd();
#else
	Point3 org = origin() + towards() * scale;
	Vector3 dx = right() * scale * tan(x_fov());
	Vector3 dy = up() * scale * tan(y_fov());
	Point3 ur = org + dx + dy;
	Point3 lr = org + dx - dy;
	Point3 ul = org - dx + dy;
	Point3 ll = org - dx - dy;

	glBegin(GL_LINE_LOOP);
	ur.load();
	ul.load();
	ll.load();
	lr.load();
	glEnd();

	glBegin(GL_LINE_LOOP);
	origin().load();
	ll.load();
	ul.load();
	glEnd();

	glBegin(GL_LINE_LOOP);
	origin().load();
	lr.load();
	ur.load();
	glEnd();

	GLfloat line_width;
	glGetFloatv(GL_LINE_WIDTH, &line_width);
	glLineWidth(2.0 * line_width);
	Point3 p1 = origin() + 1.2 * towards() * scale;
	Point3 p2 = p1 - (0.2 * dx) - (2.0 * towards());
	Point3 p3 = p1 + (0.2 * dx) - (2.0 * towards());
	glBegin(GL_LINES);
	origin().load();
	p1.load();
	p2.load();
	p1.load();
	p3.load();
	p1.load();
	glEnd();
	glLineWidth(line_width);
#endif
}


template <class T>
void tCamera<T>::init_exterior_view(const Box3 &box)
{
    //assert(!box.is_empty());
    // It's OK to be empty.  For example, a single polygon.
    value_type r = 0.5 * box.diagonal_length();
    assert(real::is_positive(r));
    Point3 o = box.centroid() + (2 * r * Vector3::VECTOR3_POSITIVE_Z);
    reset(o, Vector3::VECTOR3_NEGATIVE_Z, Vector3::VECTOR3_POSITIVE_Y,
          (T) DEFAULT_X_FOV, (T)DEFAULT_Y_FOV, (T)0.001 * r, (T)10.1 * r);
}


template <class T>
void tCamera<T>::init_interior_view(const Box3 &box)
{
    assert(!box.is_empty());
    value_type r = 0.5 * box.diagonal_length();
    assert(real::is_positive(r));
    Point3 o = box.centroid();
    reset(o, Vector3::VECTOR3_NEGATIVE_X, Vector3::VECTOR3_POSITIVE_Z,
          (T)DEFAULT_X_FOV, (T)DEFAULT_Y_FOV, (T)0.001 * r, (T)10.1 * r);
}


template <class T>
void tCamera<T>::init_image_view()
{
    _type = ORTHOGRAPHIC;
}


template <class T>
void tCamera<T>::rotate(const Line3 &l, value_type theta)
{
    Point3 o = origin();
    Vector3 v1 = towards();
    Vector3 v2 = up();
    o.rotate(l, theta);
    v1.rotate(l.direction(), theta);
    v2.rotate(l.direction(), theta);
    reset(o, v1, v2, x_fov(), y_fov(), near_distance(), far_distance());
}


template <class T>
void tCamera<T>::translate(const Vector3 &t)
{
    reset(origin() + t,
          towards(),
          up(),
          x_fov(),
          y_fov(),
          near_distance(),
          far_distance());
}


template <class T>
void tCamera<T>::move_to(const Point3 &p)
{
    reset(p,
          towards(),
          up(),
          x_fov(),
          y_fov(),
          near_distance(),
          far_distance());
}


enum { INSIDE, OUTSIDE, OVERLAP };


template <class T>
inline int plane_box_intersect(const tPlane<T> &p, const tBox3<T> &b)
{
    const tVector3<T> &n = p.normal();
    T d = p.D();
    const tPoint3<T> &b_min = b.min_point();
    const tPoint3<T> &b_max = b.max_point();
    tPoint3<T> v_min;
    tPoint3<T> v_max;
    for (unsigned i = 0; i < 3; i++) {
        if (n[i] >= 0.0) {
            v_min[i] = b_min[i];
            v_max[i] = b_max[i];
        } else {
            v_min[i] = b_max[i];
            v_max[i] = b_min[i];
        }
    }
    if (((n * v_min) + d) > 0.0) {
        return INSIDE;
    }
    if (((n * v_max) + d) >= 0.0) {
        return OVERLAP;
    }
    return OUTSIDE;
}


template <class T>
inline bool box_inside_plane(const tBox3<T> &b, const tPlane<T> &p)
{
    const tVector3<T> &n = p.normal();
    T d = p.D();
    const tPoint3<T> &b_min = b.min_point();
    const tPoint3<T> &b_max = b.max_point();
    tPoint3<T> v_min;
    for (unsigned i = 0; i < 3; i++) {
        if (n[i] >= 0.0) {
            v_min[i] = b_min[i];
        } else {
            v_min[i] = b_max[i];
        }
    }
    return ((n * v_min) + d) > 0;
}


template <class T>
inline bool box_outside_plane(const tBox3<T> &b, const tPlane<T> &p)
{
    const tVector3<T> &n = p.normal();
    T d = p.D();
    const tPoint3<T> &b_min = b.min_point();
    const tPoint3<T> &b_max = b.max_point();
    tPoint3<T> v_max;
    for (unsigned i = 0; i < 3; i++) {
        if (n[i] >= 0.0) {
            v_max[i] = b_max[i];
        } else {
            v_max[i] = b_min[i];
        }
    }
    return ((n * v_max) + d) < 0;
}


template <class T>
bool tCamera<T>::contains(const Box3 &b) const
{
    for (unsigned k = 0; k < 6; k++) {
        if (!box_inside_plane(b, _planes[k])) {
            return false;
        }
    }
    return true;
}


template <class T>
bool tCamera<T>::sees(const Box3 &b) const
{
    for (unsigned k = 0; k < 6; k++) {
        if (box_outside_plane(b, _planes[k])) {
            return false;
        }
    }
    return true;
}


template <class T>
bool tCamera<T>::sees_ignoring_near_plane(const Box3 &b) const
{
    if (box_outside_plane(b, _planes[FAR_PLANE])) {
        return false;
    }
    if (box_outside_plane(b, _planes[LEFT_PLANE])) {
        return false;
    }
    if (box_outside_plane(b, _planes[RIGHT_PLANE])) {
        return false;
    }
    if (box_outside_plane(b, _planes[TOP_PLANE])) {
        return false;
    }
    if (box_outside_plane(b, _planes[BOTTOM_PLANE])) {
        return false;
    }
    return true;
}


template <class T>
bool tCamera<T>::read(FILE *fp)
{
    float ox, oy, oz;
    float tx, ty, tz;
    float ux, uy, uz;
    float xfov, yfov;
    float n, f;
    if (fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f",
               &ox, &oy, &oz,
               &tx, &ty, &tz,
               &ux, &uy, &uz,
               &xfov, &yfov,
               &n, &f) == 13) {
        reset(Point3(ox, oy, oz),
              Vector3(tx, ty, tz),
              Vector3(ux, uy, uz),
              xfov, yfov,
              n, f);
        return true;
    } else {
        return false;
    }
}


template <class T>
void tCamera<T>::write(FILE *fp) const
{
    Point3 o = origin();
    Vector3 t = towards();
    Vector3 u = up();
    fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f\n",
            o.x(), o.y(), o.z(),
            t.x(), t.y(), t.z(),
            u.x(), u.y(), u.z(),
            x_fov(), y_fov(),
            near_distance(), far_distance());
}


template <class T>
bool tCamera<T>::read(const char *file_name)
{
    FILE *fp = xfopen(file_name, "r");
    read(fp);
    fclose(fp);
    return true;
}


template <class T>
void tCamera<T>::write(const char *file_name) const
{
    FILE *fp = xfopen(file_name, "w");
    write(fp);
    fclose(fp);
}


template <class T>
typename tCamera<T>::Point3 tCamera<T>::world_to_camera(const Point3 &p) const
{
    Point3 q(p);
    return q.affine_transform(inverse_matrix());
}


template <class T>
typename tCamera<T>::Point3 tCamera<T>::camera_to_world(const Point3 &p) const
{
    Point3 q(p);
    return q.affine_transform(matrix());
}


template <class T>
ostream &operator<<(ostream &s, const tCamera<T> &c)
{
    s << "origin: " << c.origin() << "\n";
    s << "towards: " << c.towards() << "\n";
    s << "up: " << c.up() << "\n";
    s << "x_fov:" << c.x_fov() << "\n";
    s << "y_fov:" << c.y_fov() << "\n";
    s << "near_distance:" << c.near_distance() << "\n";
    s << "far_distance:" << c.far_distance() << "\n";

#if 0
    s << "\ncorners:\n";
    for (unsigned i = 0; i < 8; i++) {
        s << "  " << i << "  " << c._corners[i] << "\n";
    }

    s << "\nplanes:\n";
    for (unsigned i = 0; i < 6; i++) {
        s << "  " << i
          << "  n: " << c._planes[i].normal()
          << "  D: " << c._planes[i].D() << "\n";
    }
#endif

    return s;
}

template ostream &operator<<(ostream &s, const tCamera<float> &c);
template ostream &operator<<(ostream &s, const tCamera<double> &c);


template <class T>
void tCamera<T>::reset_fov(value_type fovx, value_type fovy)
{
    _x_fov = fovx;
    _y_fov = fovy;
}

template class tCamera<float>;
template class tCamera<double>;

GTB_END_NAMESPACE
