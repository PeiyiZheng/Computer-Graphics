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
	GLfloat position_vec[4];
	GLfloat ambient_vec[4] = {0.0f, 0.0f, 0.0f, 1.0f};
	GLfloat diffuse_vec[4];
	GLfloat specular_vec[4];
	
	for(int i = 0; i < 3; ++i) {
		diffuse_vec[i] = specular_vec[i] = this->color[i];
		position_vec[i] = this->location[i];
	}

	diffuse_vec[3] = specular_vec[3] = position_vec[3] = 1.0f;

	glLightfv(GL_LIGHT0 + index, GL_POSITION, position_vec);
	glLightfv(GL_LIGHT0 + index, GL_AMBIENT, ambient_vec);
	glLightfv(GL_LIGHT0 + index, GL_DIFFUSE, diffuse_vec);
	glLightfv(GL_LIGHT0 + index, GL_SPECULAR, specular_vec);

	glLightf(GL_LIGHT0 + index, GL_CONSTANT_ATTENUATION, this->constAtten);
	glLightf(GL_LIGHT0 + index, GL_LINEAR_ATTENUATION, this->linearAtten);
	glLightf(GL_LIGHT0 + index, GL_QUADRATIC_ATTENUATION, this->quadAtten);

	glEnable(GL_LIGHT0 + index);
}