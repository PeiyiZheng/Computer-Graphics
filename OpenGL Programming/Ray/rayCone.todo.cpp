#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "rayScene.h"
#include "rayCone.h"


////////////////////////
//  Ray-tracing stuff //
////////////////////////
double RayCone::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
	return -1;
}

BoundingBox3D RayCone::setBoundingBox(void){
	Point3D p;
	p=Point3D(radius,height/2,radius);
	bBox=BoundingBox3D(center+p,center-p);
	return bBox;
}

//////////////////
// OpenGL stuff //
//////////////////
int RayCone::drawOpenGL(int materialIndex){
	// acceleration
	if(this->material->index != materialIndex) {
		this->material->drawOpenGL();
	}

	gluCylinder(gluNewQuadric(), this->radius, 0.0, this->height, this->openGLComplexity, this->openGLComplexity);
	GLUquadricObj *quadratic = gluNewQuadric();
	gluQuadricNormals(quadratic, GLU_FLAT);
	gluQuadricOrientation(quadratic, GLU_INSIDE);
	gluDisk(quadratic, 0, this->radius, this->openGLComplexity, this->openGLComplexity);

	return this->material->index;
}