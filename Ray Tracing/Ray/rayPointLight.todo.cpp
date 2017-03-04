#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "rayPointLight.h"
#include "rayScene.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////
Point3D RayPointLight::getDiffuse(Point3D cameraPosition,RayIntersectionInfo& iInfo){
	Point3D vec_n = iInfo.normal;
	Point3D vec_l = (iInfo.iCoordinate - this->location).unit();

	double d = vec_l.length();
	Point3D I_L = this->color / (this->constAtten + this->linearAtten * d + this->quadAtten * d * d);
	float vec_dot = vec_n.dot(vec_l);
	return vec_dot > 0 ? I_L * vec_dot : Point3D(0, 0, 0);
}

Point3D RayPointLight::getSpecular(Point3D cameraPosition,RayIntersectionInfo& iInfo){
	Point3D vec_n = iInfo.normal;
	Point3D vec_l = (iInfo.iCoordinate - this->location).unit();

	Point3D vec_v = (iInfo.iCoordinate - cameraPosition).unit();
	Point3D vec_r = (vec_n * 2.0 * vec_n.dot(vec_l) - vec_l).unit();

	double d = vec_l.length();
	Point3D I_L = this->color / (this->constAtten + this->linearAtten * d + this->quadAtten * d * d);

	float vec_dot = vec_r.dot(vec_v);
	return vec_dot > 0 ? I_L * pow(vec_dot, iInfo.material->specularFallOff) : Point3D(0, 0, 0);
}

int RayPointLight::isInShadow(RayIntersectionInfo& iInfo,RayShape* shape){
	Ray3D test_shadow;
	test_shadow.direction = (iInfo.iCoordinate - this->location).negate().unit();
	test_shadow.position = iInfo.iCoordinate + test_shadow.direction * Point3D(0.001, 0.001, 0.001); // Avoid hitting itself

	if (shape->intersect(test_shadow, iInfo, -1) > 0.0) {
		return 1;
	}

	return 0;
}

Point3D RayPointLight::transparency(RayIntersectionInfo& iInfo,RayShape* shape,Point3D cLimit){
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
void RayPointLight::drawOpenGL(int index){
}