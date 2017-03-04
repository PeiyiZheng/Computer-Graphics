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
Matrix3D getRotationMatrix(const Point3D& axis , const float& angle)
{
    double ux = axis[0], uy = axis[1], uz = axis[2];
    double cosine = cos(angle), sine = sin(angle);

    Matrix3D rotation_mat;
    rotation_mat.index(0, 0) = cosine + pow(ux, 2) * (1.0 - cosine);
    rotation_mat.index(0, 1) = uz * sine + ux * uy * (1.0 - cosine);
    rotation_mat.index(0, 2) = -uy * sine + ux * uz * (1.0 - cosine);
    rotation_mat.index(1, 0) = -uz * sine + ux * uy * (1.0 - cosine);
    rotation_mat.index(1, 1) = cosine + pow(uy, 2) * (1.0 - cosine);
    rotation_mat.index(1, 2) = ux * sine + uy * uz * (1.0 - cosine);
    rotation_mat.index(2, 0) = uy * sine + ux * uz * (1.0 - cosine);
    rotation_mat.index(2, 1) = -ux * sine + uy * uz * (1.0 - cosine);
    rotation_mat.index(2, 2) = cosine + pow(uz, 2) * (1.0 - cosine);

    return rotation_mat;
}

Matrix3D::Matrix3D(const Point3D& e){
	Matrix3D temp = getRotationMatrix(Point3D(0.0, 0.0, 1.0), e[2]);
	temp *= getRotationMatrix(Point3D(0.0, 1.0, 0.0), e[1]);
	temp *= getRotationMatrix(Point3D(1.0, 0.0, 0.0), e[0]);

	(*this) = temp;
}

Matrix3D::Matrix3D(const Quaternion& q){
	Matrix3D temp;

	double a = q.real;
	double b = q.imag[0], c = q.imag[1], d = q.imag[2];
	temp.index(0, 0) = 1.0 - 2.0 * (c * c + d * d);
	temp.index(1, 0) = 2.0 * (b * c - a * d);
	temp.index(2, 0) = 2.0 * (b * d + a * c);
	temp.index(0, 1) = 2.0 * (b * c + a * d);
	temp.index(1, 1) = 1.0 - 2.0 * (b * b + d * d);
	temp.index(2, 1) = 2.0 * (c * d - a * b);
	temp.index(0, 2) = 2.0 * (b * d - a * c);
	temp.index(1, 2) = 2.0 * (c * d + a * b);
	temp.index(2, 2) = 1.0 - 2.0 * (b * b + c * c);
	(*this) = temp;
}

Matrix3D Matrix3D::closestRotation(void) const {
	Matrix3D temp = (*this);
	Matrix3D r1, diagonal, r2;
	temp.SVD(r1, diagonal, r2);

	float flag = 1.0;
	for(int i = 0; i < 3; ++i) {
		flag *= diagonal.index(i, i);
		diagonal.index(i, i) = diagonal.index(i, i) < 0 ? -1.0 : 1.0;
	}

	if(flag < 0.0) {
		diagonal.index(2, 2) = -diagonal.index(2, 2);
	}

	temp = r1 * diagonal * r2;
	(*this) = temp;
	return (*this);
}
/* While these Exp and Log implementations are the direct implementations of the Taylor series, the Log
 * function tends to run into convergence issues so we use the other ones:*/
Matrix3D Matrix3D::Exp(const Matrix3D& m,int iter){
	
	Matrix3D sum = Matrix3D::IdentityMatrix();
	Matrix3D temp = m;
	double div = 1.0;
	for (int i = 1; i < iter; ++i) {
		sum += temp;
		temp *= m;
		div += 1.0;
		temp /= div;
	}

	return sum;
}
