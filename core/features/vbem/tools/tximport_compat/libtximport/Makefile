CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

SRCS = tximport.cpp
OBJS = $(SRCS:.cpp=.o)
LIB = libtximport.a

TEST_SRCS = test_tximport.cpp
TEST_BIN = test_tximport

.PHONY: all clean test

all: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

%.o: %.cpp tximport.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TEST_SRCS) $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_SRCS) -L. -ltximport

clean:
	rm -f $(OBJS) $(LIB) $(TEST_BIN)


