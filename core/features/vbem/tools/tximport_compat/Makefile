CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

CLI_SRC = tximport_compat.cpp
CLI_BIN = tximport_compat

LIB_DIR = libtximport
LIB = $(LIB_DIR)/libtximport.a

.PHONY: all clean test lib

all: lib $(CLI_BIN)

lib:
	$(MAKE) -C $(LIB_DIR)

$(CLI_BIN): $(CLI_SRC) $(LIB)
	$(CXX) $(CXXFLAGS) -I. -o $@ $< -L$(LIB_DIR) -ltximport

test: all
	$(MAKE) -C $(LIB_DIR) test
	@echo "=== CLI smoke test ==="
	./test_cli.sh

clean:
	$(MAKE) -C $(LIB_DIR) clean
	rm -f $(CLI_BIN)


