#include <GL/glew.h>
#include <GL/glut.h>
#include <math.h>
#include "rayCamera.h"



//////////////////
// OpenGL stuff //
//////////////////
void RayCamera::drawOpenGL(void){
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(this->position[0], this->position[1], this->position[2], this->direction[0], this->direction[1], this->direction[2], this->up[0], this->up[1], this->up[2]);
}

void RayCamera::rotateUp( Point3D center , float angle )
{
	double ux = this->up[0], uy = this->up[1], uz = this->up[2];
	double cosine = cos(angle), sine = sin(angle);
	Matrix3D rotation_mat;

	// rotation matrix from wiki
	rotation_mat.index(0, 0) = cosine + pow(ux, 2) * (1.0 - cosine);
	rotation_mat.index(0, 1) = ux * uy * (1.0 - cosine) - uz * sine;
	rotation_mat.index(0, 2) = ux * uz * (1.0 - cosine) + uy * sine;
	rotation_mat.index(1, 0) = uy * ux * (1.0 - cosine) + uz * sine;
	rotation_mat.index(1, 1) = cosine + pow(uy, 2) * (1.0 - cosine);
	rotation_mat.index(1, 2) = uy * uz * (1.0 - cosine) - ux * sine;
	rotation_mat.index(2, 0) = uz * ux * (1.0 - cosine) - uy * sine;
	rotation_mat.index(2, 1) = uz * uy * (1.0 - cosine) + ux * sine;
	rotation_mat.index(2, 2) = cosine + pow(uz, 2) * (1.0 - cosine);
    
    // change position
    this->position = this->position - center; // move to center
    this->position = rotation_mat * this->position; // rotate
    this->position = this->position + center; // move back

    // change direction -- up unchanged
    this->right = (rotation_mat * this->right).unit();
    this->direction = (rotation_mat * this->direction).unit();
}

void RayCamera::rotateRight( Point3D center , float angle )
{
    double ux = this->right[0], uy = this->right[1], uz = this->right[2];
	double cosine = cos(angle), sine = sin(angle);
	Matrix3D rotation_mat;

	// rotation matrix from wiki
	rotation_mat.index(0, 0) = cosine + pow(ux, 2) * (1.0 - cosine);
	rotation_mat.index(0, 1) = ux * uy * (1.0 - cosine) - uz * sine;
	rotation_mat.index(0, 2) = ux * uz * (1.0 - cosine) + uy * sine;
	rotation_mat.index(1, 0) = uy * ux * (1.0 - cosine) + uz * sine;
	rotation_mat.index(1, 1) = cosine + pow(uy, 2) * (1.0 - cosine);
	rotation_mat.index(1, 2) = uy * uz * (1.0 - cosine) - ux * sine;
	rotation_mat.index(2, 0) = uz * ux * (1.0 - cosine) - uy * sine;
	rotation_mat.index(2, 1) = uz * uy * (1.0 - cosine) + ux * sine;
	rotation_mat.index(2, 2) = cosine + pow(uz, 2) * (1.0 - cosine);
    
    // change position
    this->position = this->position - center; // move to center
    this->position = rotation_mat * this->position; // rotate
    this->position = this->position + center; // move back

    // change direction -- right unchanged
    this->up = (rotation_mat * this->up).unit();
    this->direction = (rotation_mat * this->direction).unit();
}

void RayCamera::moveForward(float dist){
}
void RayCamera::moveRight(float dist){
}
void RayCamera::moveUp(float dist){
}
