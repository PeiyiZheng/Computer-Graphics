#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "rayScene.h"
#include "raySphere.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////

double RaySphere::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
	double b = ray.direction.dot(ray.position - this->center) * 2.0;
	double c = pow((ray.position - this->center).length(), 2) - pow(this->radius, 2);

	double delta = b * b - 4.0 * c;

	if (delta < 0.0) {
		return -1;
	}

	double t1 = (-b + sqrt(delta)) / 2.0, t2 = (-b - sqrt(delta)) / 2.0;

	double t_closest = t1 < t2 ? t1 : t2;

	if (mx >= 0.0 && t_closest > mx) {
		return -1;
	}

	Point3D intersection = ray.position + ray.direction * t_closest;
	iInfo.iCoordinate = intersection;
	iInfo.normal = (intersection - this->center).unit();
	iInfo.material = this->material;

	// Texture mapping
	float theta = atan2(-(intersection[2] - this->center[2]), intersection[0] - this->center[0]);
	float u = (theta + PI) / (2.0 * PI);
	float gamma = acos(-(intersection[1] - this->center[1]) / this->radius);
	float v = gamma / PI;

	iInfo.texCoordinate[0] = u;
	iInfo.texCoordinate[1] = v;

	return t_closest;
}

BoundingBox3D RaySphere::setBoundingBox(void){
	Point3D rd = Point3D(this->radius, this->radius, this->radius);
	this->bBox = BoundingBox3D(this->center - rd, this->center + rd);
	return this->bBox;
}

//////////////////
// OpenGL stuff //
//////////////////
int RaySphere::drawOpenGL(int materialIndex){
	return -1;
}
