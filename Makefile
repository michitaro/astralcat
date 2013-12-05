################################################################################
# edit here
################################################################################
SFITSIO_DIR := $(HOME)/local
LIBRAW_DIR  := $(HOME)/local
GSL_DIR     := $(HOME)/local
BOOST_DIR   := $(HOME)/local/gcc47/include
CXX         := $(HOME)/local/gcc47/bin/g++
################################################################################

CXXFLAGS += -O3 -g -std=c++0x 
CXXFLAGS += -funroll-loops

# warning
CXXFLAGS += -Werror -Wall -Wno-sign-compare -Wno-parentheses

# openmp
CXXFLAGS += -fopenmp
LDFLAGS  += -fopenmp

# sfitsio
CXXFLAGS += -I$(SFITSIO_DIR)/include -DSLI__USE_CMATH
LDFLAGS  += -L$(SFITSIO_DIR)/lib -L$(SFITSIO_DIR)/lib64 -lsllib -lsfitsio

# libraw
CXXFLAGS += -I$(LIBRAW_DIR)/include
LDFLAGS  += -L$(LIBRAW_DIR)/lib -lraw

# gsl
CXXFLAGS += -I$(GSL_DIR)/include
LDFLAGS  += -L$(GSL_DIR)/lib -lgsl -lgslcblas

# boost
CXXFLAGS += -I$(BOOST_DIR)/include -DBOOST_UBLAS_SINGULAR_CHECK
LDFLAGS  += -L$(BOOST_DIR)/lib

exec := raw2fits combine isr sky stitch

all: $(exec)

astralcat.a: Region.o ds9.o SkyEstimator.o SplineSurface.o PolynomialFitter2D.o detect.o convolve.o utils.o Logger.o Source.o mosaic.o stack.o
	$(AR) rcs $@ $^

$(exec): %: %.o astralcat.a
	LD_RUN_PATH=$(HOME)/local/lib64:$(HOME)/local/lib:$(HOME)/local/gcc47/lib64 $(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp *.h
	$(CXX) -c $(CXXFLAGS) $<

clean:
	-rm -f *.o
	-rm -f $(exec) astralcat.a
