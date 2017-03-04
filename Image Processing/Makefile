ASSIGNMENT1_TARGET=Assignment1

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
	Util/cmdLineParser.cpp

IMAGE_SOURCE = \
	Image/bmp.cpp \
	Image/image.cpp \
	Image/image.todo.cpp \
	Image/jpeg.cpp \
	Image/lineSegments.cpp \
	Image/lineSegments.todo.cpp

ASSIGNMENT1_SOURCE=$(JPEG_SOURCE) $(UTIL_SOURCE) $(IMAGE_SOURCE) main.cpp

CFLAGS += -fpermissive -fopenmp -Wno-deprecated -Wno-write-strings -msse2
LFLAGS +=

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math
LFLAGS_RELEASE = -O3 

SRC = ./
BIN = Bin/Linux/
BIN_O = ./
INCLUDE = /usr/include/ -I. -IInclude -ISource

CC=gcc
CXX=g++
MD=mkdir

ASSIGNMENT1_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(ASSIGNMENT1_SOURCE))))


all: CFLAGS += $(CFLAGS_DEBUG)
all: LFLAGS += $(LFLAGS_DEBUG)
all: $(BIN)$(ASSIGNMENT1_TARGET)

release: CFLAGS += $(CFLAGS_RELEASE)
release: LFLAGS += $(LFLAGS_RELEASE)
release: $(BIN)$(ASSIGNMENT1_TARGET)

clean:
	rm -f $(BIN)$(ASSIGNMENT1_TARGET)

$(BIN):
	mkdir -p $(BIN)

$(BIN)$(ASSIGNMENT1_TARGET): $(ASSIGNMENT1_OBJECTS)
	mkdir -p $(BIN)
	$(CXX) -o $@ $(ASSIGNMENT1_OBJECTS) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

