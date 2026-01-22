#!/bin/bash

# This script calculates a "signal" to "total" ratio from a two-column input file.
# "Total" is the sum of the second column for all lines.
# "Signal" is the sum of the second column where the first column matches a given pattern.
#
# Usage:
#   ./calculate_ratio.sh <pattern1> [pattern2] ... < input_file.txt
#
# Example:
#   Given a file 'data.txt':
#   A 10
#   B 20
#   C 30
#   A 5
#
#   Command:
#   ./calculate_ratio.sh A < data.txt
#
#   Output:
#   Signal: 15
#   Total:  65
#   Ratio:  0.230769

# --- Script Start ---

# Check if at least one pattern was provided on the command line.
if [ "$#" -lt 1 ]; then
    echo "Error: No patterns provided." >&2
    echo "Usage: $0 <pattern1> [pattern2] ... < input_file.txt" >&2
    exit 1
fi

# Use awk to process the input from stdin.
# "$@" passes all the script's command-line arguments directly to awk.
awk '
# BEGIN block runs once before processing any lines.
BEGIN {
    # Store all command-line arguments (the patterns) into an associative array
    # for quick lookups. The value `1` is arbitrary; we just care about the key.
    for (i = 1; i < ARGC; i++) {
        patterns[ARGV[i]] = 1
    }

    # Reset ARGC to 1. This is a crucial step to prevent awk from
    # trying to open the patterns as input files after it finishes with stdin.
    ARGC = 1

    # Initialize our counters.
    signal = 0
    total = 0
}

# This main block runs for every line in the input.
{
    # We only process lines where the second column ($2) is a valid number.
    # This regex matches optional leading +/- sign, digits, and an optional decimal part.
    if ($2 ~ /^[+-]?[0-9]+([.][0-9]*)?$/) {

        # Add the value from the second column to the running total.
        total += $2

        # Check if the first column ($1) exists as a key in our patterns array.
        if ($1 in patterns) {
            # If it matches, add the value to the signal total as well.
            signal += $2
        }
    }
}

# END block runs once after all lines have been processed.
END {
    # Avoid a division-by-zero error if the total is 0.
    if (total == 0) {
        ratio = "N/A (total is zero)"
    } else {
        ratio = signal / total
    }

    # Print the final results in a clean, formatted way.
    printf "Signal: %s\n", signal
    printf "Total:  %s\n", total
    printf "Ratio:  %s\n", ratio
}
' "$@"
