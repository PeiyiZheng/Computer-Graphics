#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "rayScene.h"
#include "rayBox.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////
double RayBox::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
	// reference: http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
	Point3D v_min = (this->center - this->length * 0.5 - ray.position) / ray.direction;
	Point3D v_max = (this->center + this->length * 0.5 - ray.position) / ray.direction;

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

	if (i_min < 0) {
		return -1;
	}
	else {
		iInfo.material = this->material;
		iInfo.iCoordinate = ray.position + ray.direction * i_min;
		iInfo.normal = (iInfo.iCoordinate - this->center).unit();

		return i_min;
	}
	return -1;
	
}

BoundingBox3D RayBox::setBoundingBox(void){
	this->bBox = BoundingBox3D(this->center - this->length / 2, this->center + this->length / 2);
	return this->bBox;
}


//////////////////
// OpenGL stuff //
//////////////////
int RayBox::drawOpenGL(int materialIndex){
	return -1;
}
