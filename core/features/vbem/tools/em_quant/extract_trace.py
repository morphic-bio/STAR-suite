#!/usr/bin/env python3
"""
Extract specific iterations or ECs from a trace file.

Usage:
    extract_trace.py <trace_file> [--iter <N>] [--transcript <txp>] [--ec <ec_id>]
"""

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Extract specific parts of trace file')
    parser.add_argument('trace_file', help='Trace file to extract from')
    parser.add_argument('--iter', type=int, help='Extract specific iteration')
    parser.add_argument('--transcript', help='Extract specific transcript')
    parser.add_argument('--ec', type=int, help='Extract specific EC ID')
    parser.add_argument('--header', action='store_true', help='Include header lines')
    
    args = parser.parse_args()
    
    with open(args.trace_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Print header lines if requested
            if line.startswith('#') or (line.startswith('iter') and args.header):
                if args.header:
                    print(line)
                continue
            
            fields = line.split('\t')
            
            # Filter per-iteration lines
            if fields[0].isdigit() or fields[0] == '-1':
                iter_num = int(fields[0])
                transcript = fields[1] if len(fields) > 1 else ''
                
                match_iter = args.iter is None or iter_num == args.iter
                match_txp = args.transcript is None or transcript == args.transcript
                
                if match_iter and match_txp:
                    print(line)
            
            # Filter EC-level lines
            elif fields[0] == 'EC':
                iter_num = int(fields[1])
                ec_id = int(fields[2])
                transcript = fields[3] if len(fields) > 3 else ''
                
                match_iter = args.iter is None or iter_num == args.iter
                match_txp = args.transcript is None or transcript == args.transcript
                match_ec = args.ec is None or ec_id == args.ec
                
                if match_iter and match_txp and match_ec:
                    print(line)

if __name__ == '__main__':
    main()
