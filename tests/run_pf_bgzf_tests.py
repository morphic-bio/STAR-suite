#!/usr/bin/env python3
"""One execution per distinct arm. Fresh output required; no implicit retries."""
import argparse, gzip, hashlib, json, os, struct, subprocess, zlib
from pathlib import Path

def bgzf(data, chunk):
    out = bytearray()
    for off in range(0, len(data), chunk):
        part = data[off:off+chunk]
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = c.compress(part) + c.flush()
        out += bytes.fromhex('1f8b08040000000000ff060042430200')
        out += struct.pack('<H', len(payload)+25) + payload
        out += struct.pack('<II', zlib.crc32(part), len(part))
    out += bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
    return bytes(out)

def digest(records):
    value = zlib.crc32(b"".join(b"\n".join(fields)+b"\n" for row in zip(*records) for fields in row))
    return format(value, "08x")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--baseline', type=Path, required=True); ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--cli', type=Path, help='also test standalone multi-sample error propagation')
    ap.add_argument('--resume', action='store_true', help='reuse completed arm logs without executing them again')
    a=ap.parse_args(); a.out.mkdir(parents=True, exist_ok=a.resume)
    inp=a.out/'inputs';inp.mkdir(exist_ok=a.resume)
    bcs=['AAACCCAAGAAACCAT','AAACCCAAGAAACCCA','AAACCCAAGAAACCCT']
    fs=['ATCGATCGATCGATCG','GCTAGCTAGCTAGCTA','TTAATTAATTAATTAA']
    (inp/'whitelist.txt').write_text('\n'.join(bcs)+'\n')
    (inp/'features.csv').write_text('name,sequence\n'+''.join(f'Feature{i},{s}\n' for i,s in enumerate(fs)))
    records=[[],[],[]]; raw=[b'',b'',b'']
    for i in range(2049):
        umi=''.join('ACGT'[(i//2 >> (j*2)) & 3] for j in range(12))
        seqs=[bcs[i%3]+umi, (fs[i%3] if i%13 else 'C'*16)+'A'*12,fs[(i+1)%3]+'T'*5]
        for mate,seq in enumerate(seqs):
            name=f'@read{i}/{mate+1}'.encode();seq=seq.encode();qual=b'I'*len(seq)
            records[mate].append((name[1:],seq,qual))
            raw[mate]+=name+b' extra header fields\n'+seq+b'\n+\n'+qual+b'\n'
    for i in range(3):
        (inp/f'R{i+1}.gz').write_bytes(bgzf(raw[i], [211,379,113][i]))
        (inp/f'ordinary{i+1}.gz').write_bytes(gzip.compress(raw[i],mtime=0))
    (inp/'mismatch.gz').write_bytes(bgzf(raw[1].replace(b'@read70/2',b'@other70/2'),317))
    (inp/'short.gz').write_bytes(bgzf(b'\n'.join(raw[1].split(b'\n')[:-5])+b'\n',211))
    (inp/'truncated.gz').write_bytes(bgzf(raw[1][:-5],211))
    (inp/'long.gz').write_bytes(bgzf(b'@read0/2\n'+b'A'*700+b'\n+\n'+b'I'*700+b'\n',100))
    (inp/'longname.gz').write_bytes(bgzf(b'@'+b'n'*700+b'\nAA\n+\nII\n',127))
    damaged=bytearray(bgzf(raw[1],379));size=int.from_bytes(damaged[16:18],'little')+1;damaged[size-8]^=1
    (inp/'crc.gz').write_bytes(damaged)
    (inp/'noeof.gz').write_bytes(bgzf(raw[1],311)[:-28])
    (inp/'crlf.gz').write_bytes(bgzf(raw[1].replace(b'\n',b'\r\n'),179))
    manifest={str(p.name):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in inp.iterdir()}
    (a.out/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
    seen=set(); ledger=json.loads((a.out/'executions.json').read_text()) if a.resume else []
    def run(label,argv,success=True,extra_env=None):
        env=dict(os.environ, **(extra_env or {}))
        binary=Path(argv[0]); identity=(hashlib.sha256(binary.read_bytes()).hexdigest(),tuple(argv[1:]), tuple(sorted((extra_env or {}).items())))
        assert identity not in seen;seen.add(identity)
        previous=next((row for row in ledger if row['label']==label),None)
        if previous:
            assert previous['argv']==list(map(str,argv)) and previous['binary_sha256']==identity[0] and previous.get('extra_env',{})==(extra_env or {})
            assert (previous['exit']==0)==success
            return (a.out/f'{label}.stdout').read_text()
        with (a.out/f'{label}.stdout').open('w') as out, (a.out/f'{label}.stderr').open('w') as err:
            result=subprocess.run(argv,stdout=out,stderr=err,timeout=90,env=env)
        ledger.append({'label':label,'argv':list(map(str,argv)),'binary_sha256':identity[0],'exit':result.returncode,'extra_env':extra_env or {}})
        (a.out/'executions.json').write_text(json.dumps(ledger,indent=2))
        assert (result.returncode==0)==success, (label,result.returncode,(a.out/f'{label}.stderr').read_text()[-2000:])
        return (a.out/f'{label}.stdout').read_text()
    def decode(label,threads,permits,limit,paths,success=True):
        return run(label,[str(a.binary),'decode',str(threads),str(permits),str(limit)]+[str(inp/p) for p in paths],success)
    for label,threads,permits,paths in [('paired_sync',0,0,['R1.gz','R2.gz']),('paired_parallel',4,1,['R1.gz','R2.gz']),('triple',5,2,['R1.gz','R2.gz','R3.gz'])]:
        out=decode(label,threads,permits,0,paths)
        assert 'records=2049' in out and f'crc32={digest(records[:len(paths)])}' in out and 'balanced=1' in out,out
    run('worker_exception',[str(a.binary),'decode','4','1','0',str(inp/'R1.gz'),str(inp/'R2.gz')],False,{'PF_TEST_THROW_INFLATE':'1'})
    assert next(row for row in ledger if row['label']=='worker_exception')['exit']==1
    assert 'injected inflate failure' in (a.out/'worker_exception.stderr').read_text()
    decode('early_close',4,1,3,['R1.gz','R2.gz'])
    out=decode('crlf',3,1,0,['R1.gz','crlf.gz']); assert f'crc32={digest(records[:2])}' in out
    # Canonical BGZF reader permits a missing terminal EOF member; preserve that policy.
    decode('no_eof_marker',2,1,0,['R1.gz','noeof.gz'])
    for bad in ['mismatch','short','truncated','long','longname','crc']:
        decode('reject_'+bad,4,1,0,['R1.gz',bad+'.gz'],False)
    def assign(label,binary,mode,threads,permits,paths,success=True):
        out=a.out/label
        return run(label,[str(binary),'assign',mode,str(threads),str(permits),str(inp/'whitelist.txt'),str(inp/'features.csv'),str(out)]+[str(inp/p) for p in paths],success)
    assign('baseline',a.baseline,'off',0,0,['R1.gz','R2.gz'])
    assign('native',a.binary,'range',4,1,['R1.gz','R2.gz'])
    assign('disabled',a.binary,'off',0,0,['R1.gz','R2.gz'])
    assign('ordinary_auto',a.binary,'auto',4,1,['ordinary1.gz','ordinary2.gz'])
    assign('mixed_pair',a.binary,'auto',4,1,['R1.gz','ordinary2.gz'])
    assign('mixed_lanes',a.binary,'auto',4,1,['R1.gz','R2.gz','ordinary1.gz','ordinary2.gz'])
    assign('mixed_lanes_control',a.baseline,'off',0,0,['R1.gz','R2.gz','ordinary1.gz','ordinary2.gz'])
    assign('forced_ordinary',a.binary,'range',4,1,['ordinary1.gz','ordinary2.gz'],False)
    assign('assignment_crc',a.binary,'range',4,1,['R1.gz','crc.gz'],False)
    def counts(label):
        d=a.out/label/'sample'
        if not (d/'matrix.mtx').exists(): return None
        features=(d/'features.txt').read_text().splitlines()
        barcodes=(d/'barcodes.txt').read_text().splitlines()
        assert len(set(features))==len(features) and len(set(barcodes))==len(barcodes)
        lines=[x for x in (d/'matrix.mtx').read_text().splitlines() if not x.startswith('%')]
        rows,cols,nnz=map(int,lines[0].split()); assert (rows,cols)==(len(features),len(barcodes))
        entries={}
        for line in lines[1:]:
            i,j,value=map(int,line.split()); key=(features[i-1],barcodes[j-1])
            assert key not in entries; entries[key]=value
        assert len(entries)==nnz
        return set(features),set(barcodes),entries,(d/'stats.txt').read_text()
    base=counts('baseline');assert base
    for label in ['native','disabled','ordinary_auto','mixed_pair']:
        assert counts(label)==base, f'count/axis/diagnostic mismatch: {label}'
    assert counts('mixed_lanes')==counts('mixed_lanes_control')
    assert counts('assignment_crc') is None and counts('forced_ordinary') is None
    if a.cli:
        for label,feature in [('bad','crc.gz'),('good','R2.gz')]:
            d=a.out/label;d.mkdir(exist_ok=a.resume)
            for mate,source in [(1,'R1.gz'),(2,feature)]:
                link=d/f'sample_R{mate}_001.fastq.gz'
                if not link.exists():link.symlink_to((inp/source).resolve())
        dest=a.out/'multisample_output'
        run('multisample_failure',[str(a.cli),'-w',str(inp/'whitelist.txt'),'-f',str(inp/'features.csv'),'-d',str(dest),'-t','1','-c','2','-S','1','-R','12','--skip_emptydrops','--skip_qc_outputs','--readFilesBgzfMode','range','--bgzfReaderThreads','4',str(a.out/'bad'),str(a.out/'good')],False)
        assert next(row for row in ledger if row['label']=='multisample_failure')['exit']==1
        assert len(list(dest.rglob('matrix.mtx')))==1
    (a.out/'PASS.json').write_text(json.dumps({'passed':len(ledger),'exact_output_arms':5,'expected_records':2049},indent=2))
    print(f'PASS: {len(ledger)} distinct arms; exact matrices and axes; shared single-permit progress')
if __name__=='__main__':main()
