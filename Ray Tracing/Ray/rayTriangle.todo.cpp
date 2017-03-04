#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rayTriangle.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////
void RayTriangle::initialize(void){
	this->plane = Plane3D(this->v[0]->position, this->v[1]->position, this->v[2]->position);
}


double RayTriangle::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
	Point3D vec_n = this->plane.normal.unit();

	double t = -(this->plane.distance + ray.position.dot(vec_n)) / (ray.direction.dot(vec_n));
	if (t < 0.0 || (mx > 0.0 && t > mx)) {
		return -1;
	}

	Point3D p = ray.position + ray.direction * t;
	Point3D temp = p - ray.position;

	int idx[3][2] = {{0, 1}, {1, 2}, {2, 0}};
	for (int i = 0; i < 3; ++i) {
		this->v1 = this->v[idx[i][0]]->position - ray.position;
		this->v2 = this->v[idx[i][1]]->position - ray.position;

		Point3D n1 = this->v2 ^ this->v1;
		if (temp.dot(n1.unit()) < 0.0) {
			return -1;
		}
	}

	iInfo.iCoordinate = p;
	iInfo.material = this->material;
	//iInfo.normal = this->plane.normal.unit();

	// texture	
	Point3D edge1 = v[1]->position - v[0]->position;
	Point3D edge2 = v[2]->position - v[0]->position;

	// scan conversion
	double area = (edge1 ^ edge2).length();
	double gamma = ((p - v[0]->position) ^ edge1).length() / area;
	double beta = ((p - v[0]->position) ^ edge2).length() / area;
	double alpha = 1.0 - gamma - beta;

	iInfo.texCoordinate = v[2]->texCoordinate * gamma + v[1]->texCoordinate * beta + v[0]->texCoordinate * alpha;

	iInfo.normal = v[2]->normal * gamma;
	iInfo.normal += v[1]->normal * beta;
	iInfo.normal += v[0]->normal * alpha;
	iInfo.normal = iInfo.normal.unit();
	
	return t;
}

BoundingBox3D RayTriangle::setBoundingBox(void){
	Point3D vertexs[3] = {this->v[0]->position, this->v[1]->position, this->v[2]->position};
	this->bBox = BoundingBox3D(vertexs, 3);
	return bBox;
}

//////////////////
// OpenGL stuff //
//////////////////
int RayTriangle::drawOpenGL(int materialIndex){
	return -1;
}
