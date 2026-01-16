#!/bin/bash
# Test STAR-Flex cutadapt 3.2 compatibility mode against generated fixtures
#
# Usage: ./tests/test_cutadapt32_parity.sh
#
# This script:
# 1. Reads fixtures from tests/fixtures/cutadapt32_adapter_fixtures.json
# 2. Runs trimvalidate with --compat Cutadapt3 mode
# 3. Compares trim positions

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TRIMVALIDATE="$PROJECT_ROOT/tools/trimvalidate/trimvalidate"
FIXTURES_FILE="$PROJECT_ROOT/tests/fixtures/cutadapt32_adapter_fixtures.json"

if [ ! -x "$TRIMVALIDATE" ]; then
    echo "ERROR: trimvalidate not found at $TRIMVALIDATE"
    echo "Build it first: cd tools/trimvalidate && make"
    exit 1
fi

if [ ! -f "$FIXTURES_FILE" ]; then
    echo "ERROR: Fixtures file not found at $FIXTURES_FILE"
    echo "Generate it first: docker run --rm -v \$(pwd):/work -w /work python:3.9-slim bash -c 'pip install -q cutadapt==3.2 && python3 tests/fixtures/generate_cutadapt32_fixtures.py'"
    exit 1
fi

echo "Testing cutadapt 3.2 compatibility mode"
echo "========================================"
echo ""

PASSED=0
FAILED=0
TOTAL=0

# Parse JSON and run tests (using Python for JSON parsing)
FIXTURES_FILE="$FIXTURES_FILE" TRIMVALIDATE="$TRIMVALIDATE" python3 << 'EOF'
import json
import subprocess
import sys
import os

fixtures_file = os.environ.get('FIXTURES_FILE', 'tests/fixtures/cutadapt32_adapter_fixtures.json')
trimvalidate = os.environ.get('TRIMVALIDATE', 'tools/trimvalidate/trimvalidate')

with open(fixtures_file) as f:
    fixtures = json.load(f)

passed = 0
failed = 0

for fixture in fixtures:
    name = fixture['name']
    read = fixture['read']
    adapter = fixture['adapter']
    min_overlap = fixture['min_overlap']
    max_error_rate = fixture['max_error_rate']
    expected_trim_pos = fixture['expected_trim_pos']
    
    # Create a fake quality string (all 'I' = Q40)
    qual = 'I' * len(read)
    
    # Run trimvalidate in compat mode
    cmd = [
        trimvalidate,
        '--adapter', adapter,
        '--quality', '0',  # No quality trimming
        '--length', '1',   # Min length 1
        '--min-overlap', str(min_overlap),
        '--compat', 'Cutadapt3',
        '--single'
    ]
    
    try:
        result = subprocess.run(
            cmd,
            input=f"{read}\n{qual}\n",
            capture_output=True,
            text=True
        )
        
        # Parse output - look for "TrimPos: <number>"
        output_lines = result.stdout.strip().split('\n')
        
        actual_trim_pos = None
        for line in output_lines:
            if line.startswith('TrimPos: '):
                actual_trim_pos = int(line[9:])
                break
        
        if actual_trim_pos is None:
            print(f"  ERROR: {name}: Could not parse output")
            print(f"        stdout: {result.stdout}")
            print(f"        stderr: {result.stderr}")
            failed += 1
            continue
        
        if actual_trim_pos == expected_trim_pos:
            print(f"  PASS: {name}")
            print(f"        Expected trim at {expected_trim_pos}, got {actual_trim_pos}")
            passed += 1
        else:
            print(f"  FAIL: {name}")
            print(f"        Expected trim at {expected_trim_pos}, got {actual_trim_pos}")
            print(f"        Read: {read}")
            failed += 1
            
    except Exception as e:
        print(f"  ERROR: {name}: {e}")
        import traceback
        traceback.print_exc()
        failed += 1

print("")
print(f"Results: {passed} passed, {failed} failed out of {passed + failed} tests")

if failed > 0:
    sys.exit(1)
EOF

echo ""
echo "Test complete."
