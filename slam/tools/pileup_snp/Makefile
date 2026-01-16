ROOT_DIR := $(abspath ../../../..)
SLAM_SRC_DIR := $(ROOT_DIR)/slam/source
LEGACY_SRC_DIR := $(ROOT_DIR)/core/legacy/source

CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++11 -I$(SLAM_SRC_DIR) -I$(LEGACY_SRC_DIR)
LDFLAGS = -L$(LEGACY_SRC_DIR)/htslib -lhts -lz -lpthread

TARGET = pileup_snp
SOURCES = pileup_snp.cpp binomial_math.cpp
HTSLIB = $(LEGACY_SRC_DIR)/htslib/libhts.a

all: $(TARGET)

$(TARGET): $(SOURCES) $(HTSLIB)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o

install: $(TARGET)
	cp $(TARGET) $(ROOT_DIR)/bin/

.PHONY: all clean install

$(HTSLIB):
	$(MAKE) -C $(LEGACY_SRC_DIR)/htslib lib-static
