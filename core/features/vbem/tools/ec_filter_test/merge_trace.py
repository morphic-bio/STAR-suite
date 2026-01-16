#!/usr/bin/env python3
import sys

def parse_fields(payload):
    fields = {}
    for part in payload.split(';'):
        if not part or '=' not in part:
            continue
        k, v = part.split('=', 1)
        fields[k] = v
    return fields

def parse_list(s, cast=float):
    if not s:
        return []
    return [cast(x) for x in s.split(',') if x != '']

def merge_traces(in_path, out_path):
    merged = {}
    with open(in_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            qname, payload = line.split('\t', 1)
            fields = parse_fields(payload)

            txp_ids = parse_list(fields.get('txpIDs',''), int)
            log_compat = parse_list(fields.get('logCompat',''), float)
            aux_prob = parse_list(fields.get('auxProb',''), float)
            best_as = fields.get('bestAS')
            best_as = int(best_as) if best_as not in (None, '') else None

            entry = merged.setdefault(qname, {'bestAS': None, 'txp': {}})
            if best_as is not None:
                entry['bestAS'] = best_as if entry['bestAS'] is None else max(entry['bestAS'], best_as)

            # Merge per-txp fields
            for i, tid in enumerate(txp_ids):
                lc = log_compat[i] if i < len(log_compat) else None
                ap = aux_prob[i] if i < len(aux_prob) else None
                cur = entry['txp'].get(tid, {'logCompat': None, 'auxProb': None})

                if lc is not None:
                    cur['logCompat'] = lc if cur['logCompat'] is None else max(cur['logCompat'], lc)
                if ap is not None:
                    cur['auxProb'] = ap if cur['auxProb'] is None else max(cur['auxProb'], ap)

                entry['txp'][tid] = cur

    with open(out_path, 'w') as out:
        for qname in sorted(merged.keys()):
            entry = merged[qname]
            tids = sorted(entry['txp'].keys())
            logc = [entry['txp'][tid]['logCompat'] for tid in tids]
            auxp = [entry['txp'][tid]['auxProb'] for tid in tids]

            # stringify
            tids_s = ",".join(str(t) for t in tids)
            logc_s = ",".join("" if v is None else str(v) for v in logc)
            auxp_s = ",".join("" if v is None else str(v) for v in auxp)
            bestas_s = "" if entry['bestAS'] is None else str(entry['bestAS'])

            out.write(
                f"{qname}\t"
                f"txpIDs={tids_s};bestAS={bestas_s};logCompat={logc_s};auxProb={auxp_s};\n"
            )

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: merge_trace.py <input_trace.tsv> <output_trace.tsv>", file=sys.stderr)
        sys.exit(1)
    merge_traces(sys.argv[1], sys.argv[2])
