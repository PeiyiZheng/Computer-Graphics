#include <stdlib.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "rayGroup.h"
#include "rayTriangle.h"
#include "raySphere.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////

double RayGroup::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
	double current_closest = mx == -1 ? -1.0 : mx;

	Point3D original_position = ray.position;

	ray.direction = (this->getInverseMatrix().multDirection(ray.direction)).unit();
	ray.position = this->getInverseMatrix().multPosition(ray.position);

	
	int idx = 0;
	for (int i = 0; i < this->sNum; ++i) {
		double dist = this->shapes[i]->bBox.intersect(ray);
		if (dist != -1.0 && (mx == -1.0 || dist < mx)) {
			this->hits[idx].t = dist;
			this->hits[idx++].shape = this->shapes[i];
		}
	}

	qsort(this->hits, idx, sizeof(RayShapeHit), RayShapeHit::Compare);
	double min_t = -1.0;
	for (int i = 0; i < idx; ++i) {
		RayIntersectionInfo temp_intersect;
		double dist = this->hits[i].shape->intersect(ray, temp_intersect, mx);

		if (dist <= 0.0) {
			continue;
		}

		temp_intersect.normal = (this->getNormalMatrix().multDirection(temp_intersect.normal)).unit();
		temp_intersect.iCoordinate = this->getMatrix().multPosition(temp_intersect.iCoordinate);

		dist = (original_position - temp_intersect.iCoordinate).length();
		if (min_t == -1.0 || dist < min_t) {
			min_t = dist;
			iInfo.iCoordinate = temp_intersect.iCoordinate;
			iInfo.normal = temp_intersect.normal;
			iInfo.material = temp_intersect.material;
			iInfo.texCoordinate = temp_intersect.texCoordinate;
		}
	}
	
	// Before acceleration
	/*
	for (int i = 0; i < this->sNum; ++i) {
		RayIntersectionInfo temp_intersect;
		double dist = this->shapes[i]->intersect(ray, temp_intersect, mx);

		if (dist <= 0.0) {
			continue;
		}

		temp_intersect.normal = (this->getNormalMatrix().multDirection(temp_intersect.normal)).unit();
		temp_intersect.iCoordinate = this->getMatrix().multPosition(temp_intersect.iCoordinate);

		dist = (original_position - temp_intersect.iCoordinate).length();

		if (current_closest < 0 || dist < current_closest) {
			current_closest = dist;
			iInfo.iCoordinate = temp_intersect.iCoordinate;
			iInfo.normal = temp_intersect.normal;
			iInfo.material = temp_intersect.material;
		}
	}

	if (mx < 0 || current_closest < mx) {
		return current_closest;
	}
	
	return -1;*/
	return min_t;
}


BoundingBox3D RayGroup::setBoundingBox(void){
	Point3D* pList;
	BoundingBox3D tBBox;
	pList=new Point3D[sNum*2];
	for(int i=0;i<sNum;i++){
		tBBox=shapes[i]->setBoundingBox();
		pList[2*i  ]=tBBox.p[0];
		pList[2*i+1]=tBBox.p[1];
	}
	tBBox=BoundingBox3D(pList,sNum*2);

	delete[] pList;
	bBox=tBBox.transform(getMatrix());
	return bBox;
}

int StaticRayGroup::set(void){
	this->normalTransform = (this->localTransform.transpose()).invert();
	this->inverseTransform = this->localTransform.invert();
	return 1;
}
//////////////////
// OpenGL stuff //
//////////////////

int RayGroup::drawOpenGL(int materialIndex){
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	Matrix4D local_transform_2d = this->getMatrix();
	GLfloat local_transform_1d[16];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			local_transform_1d[i * 4 + j] = local_transform_2d.index(i, j);
		}
	}
	glMultMatrixf(local_transform_1d);
	for(int i = 0; i < this->sNum; ++i) {
		materialIndex = this->shapes[i]->drawOpenGL(materialIndex);
	}
	glPopMatrix();
	return materialIndex;
}

/////////////////////
// Animation Stuff //
/////////////////////
Matrix4D ParametrizedEulerAnglesAndTranslation::getMatrix(void){
	return Matrix4D(Matrix3D(this->value->eulerAngles), this->value->translate);
}
Matrix4D ParametrizedClosestRotationAndTranslation::getMatrix(void){
	return Matrix4D(this->value->rotation.closestRotation(), this->value->translate);
}
Matrix4D ParametrizedRotationLogarithmAndTranslation::getMatrix(void){
	return Matrix4D(Matrix3D::Exp(this->value->skewSymmetric, 25), this->value->translate);
}
Matrix4D ParametrizedQuaternionAndTranslation::getMatrix(void){
	return Matrix4D(Matrix3D(this->value->quaternion), this->value->translate);
}
