#!/bin/bash

# Script to count unique values in the second field of a file, skipping the first 3 lines.
# Usage: ./count_unique.sh <filename>

# --- Configuration ---
# The number of lines to skip from the beginning of the file.
LINES_TO_SKIP=3
# The field number to analyze for unique values.
FIELD_NUMBER=2

# --- Input Validation ---

# Check if a filename was provided as an argument.
if [ -z "$1" ]; then
  echo "Error: No filename provided."
  echo "Usage: $0 <filename>"
  exit 1
fi

FILENAME="$1"

# Check if the provided file exists and is readable.
if [ ! -f "$FILENAME" ] || [ ! -r "$FILENAME" ]; then
  echo "Error: File '$FILENAME' does not exist or is not readable."
  exit 1
fi

# --- Processing ---

# The core logic of the script is a single pipeline of commands:
#
# 1. `tail -n +$((LINES_TO_SKIP + 1)) "$FILENAME"`
#    - `tail` outputs the last part of files.
#    - The `-n +K` option outputs lines starting with the Kth line.
#    - We use `$((...))` for arithmetic expansion to calculate the starting line number.
#    - This command effectively skips the first `LINES_TO_SKIP` (3) lines.
#
# 2. `awk -v field="$FIELD_NUMBER" '{print $field}'`
#    - `awk` is a powerful pattern scanning and processing language.
#    - It processes the input line by line.
#    - `-v field="$FIELD_NUMBER"` passes the shell variable `FIELD_NUMBER` into awk as a variable named `field`.
#    - `'{print $field}'` prints the field specified by the `field` variable for each line.
#
# 3. `sort`
#    - This command sorts the lines it receives from `awk`. This is necessary
#      because `uniq` only removes *adjacent* duplicate lines.
#
# 4. `uniq`
#    - This filters out adjacent duplicate lines, leaving only one copy of each.
#
# 5. `wc -l`
#    - `wc` (word count) with the `-l` flag counts the number of lines.
#    - The final line count is the number of unique values.

echo "Processing file: $FILENAME"
echo "Skipping the first $LINES_TO_SKIP lines..."
echo "Counting unique values in field $FIELD_NUMBER..."

UNIQUE_COUNT=$(tail -n +$((LINES_TO_SKIP + 1)) "$FILENAME" | awk -v field="$FIELD_NUMBER" '{print $field}' | sort | uniq | wc -l)

# --- Output ---

echo "----------------------------------------"
echo "Total number of unique values found: $UNIQUE_COUNT"
echo "----------------------------------------"