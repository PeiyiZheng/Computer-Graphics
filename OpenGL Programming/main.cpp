#include <stdio.h>
#include <stdlib.h>
#include <Util/cmdLineParser.h>
#include <Ray/rayScene.h>
#include <Ray/rayWindow.h>
#include <Util/time.h>

static int DEFAULT_COMPLEXITY=10;
static int DEFAULT_RESOLUTION=500;

void ShowUsage(const char* c){
	printf("Usage %s:\n",c);
	printf("\t --in <Ray File>\n");
	printf("\t [--width <image width=%d>] [--height <image height=%d>]\n",DEFAULT_RESOLUTION,DEFAULT_RESOLUTION);
	printf("\t [--cplx <complexity=%d>]\n",DEFAULT_COMPLEXITY);
	printf("\t --out <Output Image>\n");
	printf("\t --rLimit <recursion limit> --cLimit <cut-off>\n");
}

int main(int argc,char* argv[]){
	RayScene scene;
	int cplx=DEFAULT_COMPLEXITY;
	int width=DEFAULT_RESOLUTION;
	int height=DEFAULT_RESOLUTION;

	cmdLineString In,Out;
	cmdLineInt Width,Height,Complexity,RLimit;
	cmdLineFloat CLimit;
	char* paramNames[]=			{"in",	"width",	"height",	"cplx",	"out",	"rLimit",	"cLimit"};
	cmdLineReadable* params[]=	{&In,	&Width,		&Height,	&Complexity,	&Out,	&RLimit,	&CLimit};


	cmdLineParse(argc-1,&argv[1],paramNames,7,params);
	if(!In.set){
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}
	if(Complexity.set){cplx=Complexity.value;}
	if(Width.set){width=Width.value;}
	if(Height.set){height=Height.value;}

	scene.read(In.value);
	if(Out.set && RLimit.set && CLimit.set) {
		scene.camera->aspectRatio=(float)width/(float)height;
		scene.group->setBoundingBox();

		Image32 img;
		double t=GetTime();
		scene.RayTrace(width,height,RLimit.value,CLimit.value,img);
		printf("Ray Tracing Time: %f\n",GetTime()-t);

		img.WriteImage(Out.value);
	}
	RayWindow::RayView(&scene,width,height,cplx);

	return EXIT_SUCCESS;
}
