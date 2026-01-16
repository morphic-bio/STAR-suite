CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++11 -I../../source -I../../source/htslib
LDFLAGS = -L../../source/htslib -lhts -lz -lpthread

TARGET = pileup_snp
SOURCES = pileup_snp.cpp binomial_math.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o

install: $(TARGET)
	cp $(TARGET) ../../bin/

.PHONY: all clean install
