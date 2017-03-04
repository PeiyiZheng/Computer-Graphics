#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "rayDirectionalLight.h"
#include "rayScene.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////


Point3D RayDirectionalLight::getDiffuse(Point3D cameraPosition,RayIntersectionInfo& iInfo){
	Point3D vec_l =this->direction.negate();

	float vec_dot = iInfo.normal.dot(vec_l);
	return vec_dot > 0 ? this->color * vec_dot : Point3D(0, 0, 0);
}


Point3D RayDirectionalLight::getSpecular(Point3D cameraPosition,RayIntersectionInfo& iInfo){
	Point3D vec_n = iInfo.normal;
	Point3D vec_l = this->direction.negate();

	Point3D vec_v = (cameraPosition - iInfo.iCoordinate).unit();
	Point3D vec_r = (vec_n * 2.0 * vec_n.dot(vec_l) - vec_l).unit();

	float vec_dot = vec_v.dot(vec_r);
	return vec_dot > 0 ? this->color * pow(vec_dot, iInfo.material->specularFallOff) : Point3D(0, 0, 0);
}


int RayDirectionalLight::isInShadow(RayIntersectionInfo& iInfo,RayShape* shape){
	Ray3D test_shadow;
	test_shadow.direction = this->direction.negate();
	test_shadow.position = iInfo.iCoordinate + test_shadow.direction * Point3D(0.001, 0.001, 0.001); // Avoid hitting itself

	if (shape->intersect(test_shadow, iInfo, -1.0) > 0.0) {
		return 1;
	}

	return 0;
}

Point3D RayDirectionalLight::transparency(RayIntersectionInfo& iInfo,RayShape* shape,Point3D cLimit){
	if (iInfo.material->transparent[0] > cLimit[0] && iInfo.material->transparent[1] > cLimit[1] && iInfo.material->transparent[2] > cLimit[2]) {
		RayIntersectionInfo temp = iInfo; // avoid iInfo been modified by intersect function
		Point3D trans = Point3D(1, 1, 1);

		return isInShadow(temp, shape) ? iInfo.material->transparent * transparency(iInfo, shape, cLimit / iInfo.material->transparent) : trans;
	}
	else {
		return Point3D(1, 1, 1);
	}
}

//////////////////
// OpenGL stuff //
//////////////////
void RayDirectionalLight::drawOpenGL(int index){
}