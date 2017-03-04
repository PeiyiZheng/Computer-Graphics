#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "rayScene.h"
#include "rayCylinder.h"


////////////////////////
//  Ray-tracing stuff //
////////////////////////
double RayCylinder::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
}

BoundingBox3D RayCylinder::setBoundingBox(void){
	Point3D p;
	p=Point3D(radius,height/2,radius);
	bBox=BoundingBox3D(center+p,center-p);
	return bBox;
}


//////////////////
// OpenGL stuff //
//////////////////
int RayCylinder::drawOpenGL(int materialIndex){
	// acceleration
	if(this->material->index != materialIndex) {
		this->material->drawOpenGL();
	}

	glPushMatrix();
	glTranslatef(this->center[0], this->center[1], this->center[2]);
	glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
	glTranslatef(0.0, 0.0, -0.5 * this->height);
	gluCylinder(gluNewQuadric(), this->radius, this->radius, this->height, this->openGLComplexity, this->openGLComplexity);

	GLUquadricObj *disk1 = gluNewQuadric(), *disk2 = gluNewQuadric();

	gluQuadricNormals(disk1, GLU_FLAT);
	gluQuadricOrientation(disk1, GLU_INSIDE);
	gluDisk(disk1, 0.0, this->radius, this->openGLComplexity, this->openGLComplexity);
	
	glTranslatef(0.0, 0.0, this->height);

	gluQuadricNormals(disk2, GLU_FLAT);
	gluQuadricOrientation(disk2, GLU_OUTSIDE);
	gluDisk(disk2, 0.0, this->radius, this->openGLComplexity, this->openGLComplexity);
	glTranslatef(-this->center[0], -this->center[1], -this->center[2]);
	glPopMatrix();
	return this->material->index;
}