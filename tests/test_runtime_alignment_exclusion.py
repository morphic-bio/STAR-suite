#!/usr/bin/env python3
"""Native integration tests for runtime BED exclusion using one immutable index.

Run serially with other native/benchmark jobs. Output directory must be fresh.
Requires Python 3 and STAR; uses synthetic sequence only, no external datasets.
"""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import random
import subprocess


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--star', required=True)
    p.add_argument('--baseline-star')
    p.add_argument('--out', required=True)
    args=p.parse_args()
    out=Path(args.out); out.mkdir(parents=True,exist_ok=False)
    rng=random.Random(4275)
    seq=lambda n: ''.join(rng.choices('ACGT',k=n))
    reverse=lambda s: s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    genome={name:seq(4000) for name in ['equalGood','weakGood','onlyMasked','edge','splice']}
    genome['equalBad']=genome['equalGood']
    genome['weakBad']=genome['weakGood']
    weak=list(genome['weakGood'])
    for pos in [540,580,620]: weak[pos]='ACGT'[('ACGT'.index(weak[pos])+1)%4]
    genome['weakGood']=''.join(weak)
    repeated=seq(4000)
    for i in range(12): genome[f'repeat{i}']=repeated
    sp=list(genome['splice']);sp[1100:1102]='GT';sp[2098:2100]='AG';genome['splice']=''.join(sp)
    (out/'genome.fa').write_text(''.join(f'>{k}\n{s}\n' for k,s in genome.items()))
    (out/'junctions.tab').write_text(''.join(f'splice\t1101\t{end}\t1\n' for end in [2100,2500,3000,3500]))
    reads={
        'equal':genome['equalGood'][500:650],
        'weak':genome['weakBad'][500:650],
        'only':genome['onlyMasked'][500:650],
        'repeat':repeated[500:650],
        'edgeLeft':genome['edge'][850:1000],
        'edgeRight':genome['edge'][1200:1350],
        'edgeCross':genome['edge'][950:1100],
        'splice':genome['splice'][1025:1100]+genome['splice'][2100:2175],
        'junctionFlank':genome['splice'][1000:1100],
    }
    reads={name+strand:(s if strand=='_f' else reverse(s)) for name,s in reads.items() for strand in ['_f','_r']}
    (out/'reads.fastq').write_text(''.join(f'@{n}\n{s}\n+\n'+('I'*len(s))+'\n' for n,s in reads.items()))
    masks=[('equalBad',0,4000),('weakBad',0,4000),('onlyMasked',0,4000),('edge',1000,1100),('edge',1080,1200),('splice',1500,1600)]
    masks += [(f'repeat{i}',0,4000) for i in range(1,12)]
    rng.shuffle(masks)
    (out/'mask.bed').write_text('browser position edge:1000-1200\ntrack name=test\n# zero based\n'+''.join(f'{c}\t{s}\t{e}\n' for c,s,e in masks))
    (out/'exon.bed').write_text((out/'mask.bed').read_text()+'splice\t1000\t1100\n')
    (out/'neutral.bed').write_text('edge\t3900\t3910\n')

    def run(name,argv,expect_ok=True):
        (out/(name+'.argv.json')).write_text(json.dumps(argv,indent=2)+'\n')
        with (out/(name+'.log')).open('w') as log: result=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT)
        (out/(name+'.exit.json')).write_text(json.dumps({'exit_status':result.returncode})+'\n')
        assert (result.returncode==0)==expect_ok,(name,result.returncode,out/(name+'.log'))

    index=out/'index';index.mkdir()
    run('index',[args.star,'--runMode','genomeGenerate','--genomeDir',str(index),'--genomeFastaFiles',str(out/'genome.fa'),
                 '--genomeSAindexNbases','6','--genomeChrBinNbits','12','--sjdbFileChrStartEnd',str(out/'junctions.tab'),
                 '--sjdbOverhang','149','--runThreadN','1','--outFileNamePrefix',str(out/'index_build.')])
    def index_hashes():
        return {f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in index.iterdir() if f.is_file()}
    before=index_hashes()
    def align(name,extra=(),binary=None,expect_ok=True):
        target=out/name;target.mkdir()
        argv=[binary or args.star,'--genomeDir',str(index),'--readFilesIn',str(out/'reads.fastq'),'--runThreadN','1',
              '--outFileNamePrefix',str(target)+'/', '--outSAMtype','SAM','--outSAMunmapped','Within',
              '--outSAMattributes','NH','HI','AS','nM','--outFilterScoreMin','100','--outFilterMatchNmin','100',
              '--outFilterScoreMinOverLread','0','--outFilterMatchNminOverLread','0',
              '--seedMultimapNmax','5','--winAnchorMultimapNmax','3','--outFilterMultimapNmax','3']+list(extra)
        run(name,argv,expect_ok)
        if not expect_ok:
            assert 'alignExcludeRegions' in (out/(name+'.log')).read_text(),name
            return
        assert (target/'Log.final.out').exists()
        records=collections.defaultdict(list)
        for line in (target/'Aligned.out.sam').read_text().splitlines():
            if not line.startswith('@'): records[line.split('\t')[0]].append(line.split('\t'))
        assert set(records)==set(reads)
        return {name:sorted(rows) for name,rows in records.items()}

    off=align('off')
    explicit=align('explicit_off',['--alignExcludeRegions','-'])
    assert off==explicit,'default/explicit off differ'
    if args.baseline_star: assert off==align('old_binary',binary=args.baseline_star),'mask-off changed native alignment records'
    neutral=align('neutral',['--alignExcludeRegions',str(out/'neutral.bed')])
    assert off==neutral,'unrelated mask changed seed multiplicity or alignment records'
    on=align('on',['--alignExcludeRegions',str(out/'mask.bed')])
    exon=align('exon',['--alignExcludeRegions',str(out/'exon.bed')])
    for strand in ['_f','_r']:
        def mapped(records,name): return [r for r in records[name+strand] if not int(r[1])&4]
        assert len(mapped(off,'equal'))==2
        assert {r[2] for r in mapped(off,'weak')}=={'weakBad'}
        assert not mapped(off,'repeat'),'repeat limit fixture did not exceed limit'
        assert not mapped(off,'junctionFlank'),'redundant SJ flank fixture did not exceed anchor limit'
        assert off['junctionFlank'+strand]==on['junctionFlank'+strand]
        for name,chrom in [('equal','equalGood'),('weak','weakGood'),('repeat','repeat0'),('edgeLeft','edge'),('edgeRight','edge'),('splice','splice')]:
            result=mapped(on,name)
            assert len(result)==1 and result[0][2]==chrom,(name,strand,result)
            assert bool(int(result[0][1])&16)==(strand=='_r')
            assert 'NH:i:1' in result[0]
        assert 'N' in mapped(on,'splice')[0][5]
        assert not mapped(on,'only')
        assert not mapped(on,'edgeCross'),'extension across mask escaped final block check'
        assert not mapped(exon,'splice'),'annotated junction seed escaped exon mask'
        for name in ['edgeLeft','edgeRight','splice']: assert off[name+strand]==on[name+strand]
    for name,content in {
        'unknown':'missing\t0\t10\n','negative':'edge\t-1\t20\n','empty_interval':'edge\t20\t20\n',
        'reversed':'edge\t30\t20\n','outside':'edge\t0\t4001\n','overflow':'edge\t0\t184467440737095516160\n',
        'decimal':'edge\t0\t1.5\n','empty_file':'# no intervals\n',
    }.items():
        bed=out/(name+'.bed');bed.write_text(content)
        align(name,['--alignExcludeRegions',str(bed)],expect_ok=False)
    align('missing_file',['--alignExcludeRegions',str(out/'absent.bed')],expect_ok=False)
    align('chimeric_guard',['--alignExcludeRegions',str(out/'mask.bed'),'--chimSegmentMin','12'],expect_ok=False)
    align('flex_guard',['--alignExcludeRegions',str(out/'mask.bed'),'--flex','yes'],expect_ok=False)
    assert before==index_hashes(),'runtime exclusion changed index files'
    (out/'COMPLETE.json').write_text(json.dumps({'status':'PASS','reads':len(reads),'index_sha256':before,
        'star_sha256':hashlib.sha256(Path(args.star).read_bytes()).hexdigest(),'random_seed':4275},indent=2)+'\n')
    print('PASS:',out)


if __name__=='__main__': main()
