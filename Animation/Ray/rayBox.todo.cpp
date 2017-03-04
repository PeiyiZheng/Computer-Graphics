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
	bBox=BoundingBox3D(center-(length/2),center+(length/2));
	return bBox;
}



//////////////////
// OpenGL stuff //
//////////////////

int RayBox::drawOpenGL(int materialIndex){
	// acceleration
	if(this->material->index != materialIndex) {
		this->material->drawOpenGL();
	}

	// codes modified from http://nehe.gamedev.net/tutorial/quadrics/20001/
	glBegin(GL_QUADS);          // Start Drawing Quads

	float bottom_c_l[3], top_c_l[3];
	for(int i = 0; i < 3; ++i) {
		bottom_c_l[i] = this->center[i] - this->length[i] * 0.5f;
		top_c_l[i] = this->center[i] + this->length[i] * 0.5f;
	}

	// Left Face
    glNormal3f(-1.0f, 0.0f, 0.0f);      // Normal Facing Left
    glVertex3f(bottom_c_l[0], bottom_c_l[1], bottom_c_l[2]);  // Bottom Left Of The Texture and Quad
    glVertex3f(bottom_c_l[0], bottom_c_l[1], top_c_l[2]);  // Bottom Right Of The Texture and Quad
    glVertex3f(bottom_c_l[0], top_c_l[1], top_c_l[2]);  // Top Right Of The Texture and Quad
    glVertex3f(bottom_c_l[0], top_c_l[1], bottom_c_l[2]);  // Top Left Of The Texture and Quad

    // Top Face
    glBegin(GL_QUADS);
    glNormal3f(0.0f, 1.0f, 0.0f);      // Normal Facing Up
    glVertex3f(bottom_c_l[0], top_c_l[1], bottom_c_l[2]);  // Top Left Of The Texture and Quad
    glVertex3f(bottom_c_l[0], top_c_l[1], top_c_l[2]);  // Bottom Left Of The Texture and Quad
    glVertex3f(top_c_l[0], top_c_l[1], top_c_l[2]);  // Bottom Right Of The Texture and Quad
    glVertex3f(top_c_l[0], top_c_l[1], bottom_c_l[2]);  // Top Right Of The Texture and Quad

    // Right face
    glNormal3f(1.0f, 0.0f, 0.0f);      // Normal Facing Right
    glVertex3f(top_c_l[0], top_c_l[1], bottom_c_l[2]);  // Bottom Right Of The Texture and Quad
    glVertex3f(top_c_l[0], top_c_l[1], top_c_l[2]);  // Top Right Of The Texture and Quad
    glVertex3f(top_c_l[0], bottom_c_l[1], top_c_l[2]);  // Top Left Of The Texture and Quad
    glVertex3f(top_c_l[0], bottom_c_l[1], bottom_c_l[2]);  // Bottom Left Of The Texture and Quad

    // Bottom Face
    glNormal3f(0.0f,-1.0f, 0.0f);      // Normal Facing Down
    glVertex3f(top_c_l[0], bottom_c_l[1], bottom_c_l[2]);  // Top Right Of The Texture and Quad
    glVertex3f(top_c_l[0], bottom_c_l[1], top_c_l[2]);  // Top Left Of The Texture and Quad
    glVertex3f(bottom_c_l[0], bottom_c_l[1], top_c_l[2]);  // Bottom Left Of The Texture and Quad
    glVertex3f(bottom_c_l[0], bottom_c_l[1], bottom_c_l[2]);  // Bottom Right Of The Texture and Quad

    // Front Face
    glNormal3f(0.0f, 0.0f, 1.0f);      // Normal Facing Forward
    glVertex3f(top_c_l[0], bottom_c_l[1], top_c_l[2]);  // Bottom Left Of The Texture and Quad
    glVertex3f(top_c_l[0], top_c_l[1], top_c_l[2]);  // Bottom Right Of The Texture and Quad
    glVertex3f(bottom_c_l[0], top_c_l[1], top_c_l[2]);  // Top Right Of The Texture and Quad
    glVertex3f(bottom_c_l[0], bottom_c_l[1], top_c_l[2]);  // Top Left Of The Texture and Quad

    // Back Face
    glNormal3f(0.0f, 0.0f,-1.0f);      // Normal Facing Away
    glVertex3f(top_c_l[0], top_c_l[1], bottom_c_l[2]);  // Bottom Right Of The Texture and Quad
    glVertex3f(top_c_l[0], bottom_c_l[1], bottom_c_l[2]);  // Top Right Of The Texture and Quad
    glVertex3f(bottom_c_l[0], bottom_c_l[1], bottom_c_l[2]);  // Top Left Of The Texture and Quad
    glVertex3f(bottom_c_l[0], top_c_l[1], bottom_c_l[2]);  // Bottom Left Of The Texture and Quad
    glEnd(); 

	return this->material->index;
}