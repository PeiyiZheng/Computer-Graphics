#include "rayScene.h"
#include <GL/glew.h>
#include <GL/glut.h>
#include <math.h>

///////////////////////
// Ray-tracing stuff //
///////////////////////

Point3D RayScene::Reflect(Point3D v,Point3D n){
	// similar to diffuse
	return (n * 2.0 * n.dot(v) - v).unit();
}

int RayScene::Refract(Point3D v, Point3D n, double ir, Point3D &refract) {
	double temp = v.dot(n);
	refract = (v - n * temp) * ir - n * sqrt(1.0 - ir * ir * (1.0 - temp * temp));

	// why do we need to return a value? any situation will cause refract fails?

	return 1;
}

Ray3D RayScene::GetRay(RayCamera* camera,int i,int j,int width,int height){
	Point3D p0 = camera->position;
	Point3D toward = camera->direction;
	Point3D up = camera->up;
	Point3D right = camera->right;

	double theta = tan(camera->heightAngle / 2.0);

	Point3D vec_h, vec_w, p1, p2;

	p1 = toward - up * theta;
	p2 = toward + up * theta;
	vec_h = (p1 - p2) * ((height - j + 0.5) / height);

	theta *= camera->aspectRatio;
	p1 = toward - right * theta;
	p2 = toward + right * theta;
	vec_w = (p2 - p1) * ((i + 0.5) / width);

	Point3D result = p1 + up * (theta / camera->aspectRatio) + vec_h + vec_w;
	Ray3D ray;
	ray.position = p0;
	ray.direction = result.unit();

	return ray;
}

Point3D RayScene::GetColor(Ray3D ray,int rDepth,Point3D cLimit){
	RayIntersectionInfo iInfo;
	double t = this->group->intersect(ray, iInfo, -1.0);

	if (t > 0.0) {
		Point3D p = Point3D(0, 0, 0);
		p = iInfo.material->emissive + iInfo.material->ambient * this->ambient;
		for (int i = 0; i < this->lightNum; ++i) {
			Point3D diffuse_use;
			if (iInfo.material->tex == NULL) {
				diffuse_use = this->lights[i]->getDiffuse(ray.position, iInfo);
			} else {
				Pixel32 temp = iInfo.material->tex->img->BilinearSample(iInfo.texCoordinate[0], iInfo.texCoordinate[1]);
				diffuse_use = Point3D(temp.r, temp.g, temp.b);
			}

			Point3D  diffuse_val = iInfo.material->diffuse * diffuse_use;
			Point3D specular_val = iInfo.material->specular * this->lights[i]->getSpecular(ray.position, iInfo);

			RayIntersectionInfo temp = iInfo;
			Point3D trans = Point3D(1, 1, 1);
			if (this->lights[i]->isInShadow(temp, this->group)) {
				trans = this->lights[i]->transparency(temp, this->group, cLimit).unit();
			}
			
			p += (diffuse_val + specular_val) * trans;
		}
		Ray3D reflect, refract;

		if (rDepth > 0) {
			if (iInfo.material->specular[0] > cLimit[0] && iInfo.material->specular[1] > cLimit[1] && iInfo.material->specular[2] > cLimit[2]) {
				reflect.direction = Reflect(ray.direction.negate(), iInfo.normal).unit();
				reflect.position = iInfo.iCoordinate + reflect.direction * Point3D(0.001, 0.001, 0.001);
				Point3D reflect_color = GetColor(reflect, rDepth - 1, cLimit / iInfo.material->specular);
				if (reflect_color[0] != background[0] && reflect_color[1] != background[1] && reflect_color[2] != background[2]) {
					p += reflect_color * iInfo.material->specular;
				}
			}

			if (iInfo.material->transparent[0] > cLimit[0] && iInfo.material->transparent[1] > cLimit[1] && iInfo.material->transparent[2] > cLimit[2]) {
				Point3D refractDir;

				double etai = 1.0, etat = iInfo.material->refind;
				if (iInfo.normal.dot(ray.direction) < 0) {
					double mid = etai;
					etai = etat;
					etat = mid;
				}
				int refract_result = Refract(ray.direction, iInfo.normal, etai / etat, refractDir);
				refract.direction = refractDir.unit();
				refract.position = iInfo.iCoordinate + refract.direction * Point3D(0.001, 0.001, 0.001);

				Point3D refract_color = GetColor(refract, rDepth - 1, cLimit / iInfo.material->transparent);
				p += refract_color * iInfo.material->transparent;
			}
		}
		
		for (int i = 0; i < 3; ++i) {
			p[i] = p[i] < 0.0 ? 0.0 : p[i];
			p[i] = p[i] > 1.0 ? 1.0 : p[i];
		}

		return p;
	}
	else {
		return this->background;
	}
}

//////////////////
// OpenGL stuff //
//////////////////
void RayMaterial::drawOpenGL(void){
	GLfloat ambient_vec[4] = {this->ambient[0], this->ambient[1], this->ambient[2], 1.0f};
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient_vec);

	GLfloat diffuse_vec[4] = {this->diffuse[0], this->diffuse[1], this->diffuse[2], 1.0f};
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_vec);

	GLfloat specular_vec[4] = {this->specular[0], this->specular[1], this->specular[2], 1.0f};
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_vec);

	GLfloat emissive_vec[4] = {this->emissive[0], this->emissive[1], this->emissive[2], 1.0f};
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emissive_vec);

	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, this->specularFallOff);

	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR);
	//glColor3f(this->transparent[0], this->transparent[1], this->transparent[2]); 

	if(this->tex != NULL) {
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, this->tex->openGLHandle);
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	}
}

void RayTexture::setUpOpenGL(void){
	unsigned char *texture = new unsigned char[this->img->width() * this->img->height() * 3];
	for(int i = 0; i < this->img->height(); ++i) {
		for(int j = 0; j < this->img->width(); ++j) {
			Pixel32& temp = this->img->pixel(j, i);
			int idx = 3 * (i * this->img->width() + j);
			texture[idx + 0] = temp.r;
			texture[idx + 1] = temp.g;
			texture[idx + 2] = temp.b;
		} 
	}

	glGenTextures(1, &this->openGLHandle);
	glBindTexture(GL_TEXTURE_2D, this->openGLHandle);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->img->width(), this->img->height(), 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&texture[0]);
	glBindTexture(GL_TEXTURE_2D, 0);
}