#include <stdlib.h>
#include <math.h>

#include <SVD/SVDFit.h>
#include <SVD/MatrixMNTC.h>

#include "geometry.h"


///////////////////////
// Ray-tracing stuff //
///////////////////////

double BoundingBox3D::intersect(const Ray3D &ray) const {

	// reference: http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
	Point3D v_min = (this->p[0] - ray.position) / ray.direction;
	Point3D v_max = (this->p[1] - ray.position) / ray.direction;

	double i_min, i_max, t_min, t_max;
	i_min = v_min[0] < v_max[0] ? v_min[0] : v_max[0];
	i_max = v_min[0] > v_max[0] ? v_min[0] : v_max[0];
	t_min = v_min[1] < v_max[1] ? v_min[1] : v_max[1];
	t_max = v_min[1] > v_max[1] ? v_min[1] : v_max[1];

	if (i_min > t_max || i_max < t_min) {
		return -1;
	}

	if (i_max < 0 || t_max < 0) {
		return -1;
	}

	i_min = i_min > t_min ? i_min : t_min;
	i_max = i_max < t_max ? i_max : t_max;
	t_min = v_min[2] < v_max[2] ? v_min[2] : v_max[2];
	t_max = v_min[2] > v_max[2] ? v_min[2] : v_max[2];

	if (t_max < 0 || t_min > i_max || i_min > t_max) {
		return -1;
	}

	i_min = i_min < t_min ? i_min : t_min;

	return i_min < 0 ? 0 : i_min;
}

/////////////////////
// Animation stuff //
/////////////////////
Matrix3D::Matrix3D(const Point3D& e){
	(*this)=Matrix3D();
}

Matrix3D::Matrix3D(const Quaternion& q){
	(*this)=Matrix3D();
}
Matrix3D Matrix3D::closestRotation(void) const {
	return (*this);
}
/* While these Exp and Log implementations are the direct implementations of the Taylor series, the Log
 * function tends to run into convergence issues so we use the other ones:*/
Matrix3D Matrix3D::Exp(const Matrix3D& m,int iter){
	return m;
}
