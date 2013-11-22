################################################################################
# edit here
################################################################################
SFITSIO_DIR := $(HOME)/local
LIBRAW_DIR  := $(HOME)/local/packages/libraw/0.15.4
CXX = /usr/local/Cellar/gcc47/4.7.3/bin/g++-4.7


CXXFLAGS += -O3 -g -std=c++0x

# sfitsio
CXXFLAGS += -I$(SFITSIO_DIR)/include
LDFLAGS  += -L$(SFITSIO_DIR)/lib -L$(SFITSIO_DIR)/lib64 -lsllib -lsfitsio

# libraw
CXXFLAGS += -I$(LIBRAW_DIR)/include
LDFLAGS  += -L$(LIBRAW_DIR)/lib -lraw

exec := raw2fits

$(exec): %: %.o
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $<
