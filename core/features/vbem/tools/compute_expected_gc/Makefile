CC = gcc
CFLAGS = -O3 -Wall -Wextra
LDFLAGS = -lm

TARGET = compute_expected_gc
SOURCE = compute_expected_gc.c

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE) $(LDFLAGS)

clean:
	rm -f $(TARGET)

install: $(TARGET)
	cp $(TARGET) ../../bin/

.PHONY: all clean install
