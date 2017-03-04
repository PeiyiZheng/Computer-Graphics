ASSIGNMENT4_TARGET=Assignment4

JPEG_SOURCE = \
	JPEG/jcapimin.cpp \
	JPEG/jcapistd.cpp \
	JPEG/jccoefct.cpp \
	JPEG/jccolor.cpp \
	JPEG/jcdctmgr.cpp \
	JPEG/jchuff.cpp \
	JPEG/jcinit.cpp \
	JPEG/jcmainct.cpp \
	JPEG/jcmarker.cpp \
	JPEG/jcmaster.cpp \
	JPEG/jcomapi.cpp \
	JPEG/jcparam.cpp \
	JPEG/jcphuff.cpp \
	JPEG/jcprepct.cpp \
	JPEG/jcsample.cpp \
	JPEG/jctrans.cpp \
	JPEG/jdapimin.cpp \
	JPEG/jdapistd.cpp \
	JPEG/jdatadst.cpp \
	JPEG/jdatasrc.cpp \
	JPEG/jdcoefct.cpp \
	JPEG/jdcolor.cpp \
	JPEG/jddctmgr.cpp \
	JPEG/jdhuff.cpp \
	JPEG/jdinput.cpp \
	JPEG/jdmainct.cpp \
	JPEG/jdmarker.cpp \
	JPEG/jdmaster.cpp \
	JPEG/jdmerge.cpp \
	JPEG/jdphuff.cpp \
	JPEG/jdpostct.cpp \
	JPEG/jdsample.cpp \
	JPEG/jdtrans.cpp \
	JPEG/jerror.cpp \
	JPEG/jfdctflt.cpp \
	JPEG/jfdctfst.cpp \
	JPEG/jfdctint.cpp \
	JPEG/jidctflt.cpp \
	JPEG/jidctfst.cpp \
	JPEG/jidctint.cpp \
	JPEG/jidctred.cpp \
	JPEG/jmemmgr.cpp \
	JPEG/jmemnobs.cpp \
	JPEG/jquant1.cpp \
	JPEG/jquant2.cpp \
	JPEG/jutils.cpp 

UTIL_SOURCE = \
	Util/cmdLineParser.cpp \
	Util/geometry.cpp \
	Util/geometry.todo.cpp \
	Util/time.cpp

RAY_SOURCE = \
	Ray/mouse.cpp \
	Ray/rayBox.cpp \
	Ray/rayBox.todo.cpp \
	Ray/rayCamera.cpp \
	Ray/rayCamera.todo.cpp \
	Ray/rayCone.cpp \
	Ray/rayCone.todo.cpp \
	Ray/rayCylinder.cpp \
	Ray/rayCylinder.todo.cpp \
	Ray/rayDirectionalLight.cpp \
	Ray/rayDirectionalLight.todo.cpp \
	Ray/rayFileInstance.cpp \
	Ray/rayGroup.cpp \
	Ray/rayGroup.todo.cpp \
	Ray/rayKey.cpp \
	Ray/rayPointLight.cpp \
	Ray/rayPointLight.todo.cpp \
	Ray/rayScene.cpp \
	Ray/rayScene.todo.cpp \
	Ray/raySphere.cpp \
	Ray/raySphere.todo.cpp \
	Ray/raySpotLight.cpp \
	Ray/raySpotLight.todo.cpp \
	Ray/rayTriangle.cpp \
	Ray/rayTriangle.todo.cpp \
	Ray/rayWindow.cpp

IMAGE_SOURCE = \
	Image/bmp.cpp \
	Image/image.cpp \
	Image/image.todo.cpp \
	Image/jpeg.cpp \
	Image/lineSegments.cpp \
	Image/lineSegments.todo.cpp

GL_SOURCE = \
	GL/glew.c

ASSIGNMENT4_SOURCE=$(JPEG_SOURCE) $(UTIL_SOURCE) $(IMAGE_SOURCE) $(RAY_SOURCE) $(GL_SOURCE) main.cpp

OS := $(shell uname -s)
ifeq ($(OS),Darwin)
CFLAGS += -fpermissive -Wno-deprecated -Wno-write-strings -msse2
LFLAGS += -framework OpenGL -framework glut
INCLUDE = /usr/include/ -I. -IInclude -ISource -I/usr/include/malloc
else
CFLAGS += -fpermissive -fopenmp -Wno-deprecated -Wno-write-strings -msse2
LFLAGS += -lGL -lGLU -lglut
INCLUDE = /usr/include/ -I. -IInclude -ISource
endif

CFLAGS_DEBUG = -DDEBUG -g3 -DGLEW_STATIC
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math
LFLAGS_RELEASE = -O3 

SRC = ./
BIN = Bin/Linux/
BIN_O = ./

CC=gcc
CXX=g++
MD=mkdir

ASSIGNMENT4_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(ASSIGNMENT4_SOURCE))))


all: CFLAGS += $(CFLAGS_DEBUG)
all: LFLAGS += $(LFLAGS_DEBUG)
all: $(BIN)$(ASSIGNMENT4_TARGET)

release: CFLAGS += $(CFLAGS_RELEASE)
release: LFLAGS += $(LFLAGS_RELEASE)
release: $(BIN)$(ASSIGNMENT4_TARGET)

clean:
	rm -f $(BIN)$(ASSIGNMENT4_TARGET)
	rm -f GL/*.o Image/*.o Ray/*.o SVD/*.o Util/*.o main.o

$(BIN):
	mkdir -p $(BIN)

$(BIN)$(ASSIGNMENT4_TARGET): $(ASSIGNMENT4_OBJECTS)
	mkdir -p $(BIN)
	$(CXX) -o $@ $(ASSIGNMENT4_OBJECTS) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

