
import sys
import mappy as mp
import edlib
from .utils import zopen, check_inputs


degenNuc = [("R", "A"), ("R", "G"),
            ("M", "A"), ("M", "C"),
            ("W", "A"), ("W", "T"),
            ("S", "C"), ("S", "G"),
            ("Y", "C"), ("Y", "T"),
            ("K", "G"), ("K", "T"),
            ("V", "A"), ("V", "C"), ("V", "G"),
            ("H", "A"), ("H", "C"), ("H", "T"),
            ("D", "A"), ("D", "G"), ("D", "T"),
            ("B", "C"), ("B", "G"), ("B", "T"),
            ("N", "G"), ("N", "A"), ("N", "T"), ("N", "C"),
            ("X", "G"), ("X", "A"), ("X", "T"), ("X", "C")]


def find_all_splice_GT(seq, offset=6):
    # check reference end for splice donor, might be off by a few bp when small exons missed
    possible = []
    for i in range(0, offset):
        pos = offset-i
        if seq[pos:pos+2].startswith('GT'):
            possible.append(i)
    return possible


def find_all_splice_AG(seq, offset=6):
    # check reference end for splice acceptor, might be off by a few bp when small exons missed
    possible = []
    for i in range(0, offset):
        pos = offset-i
        if seq[:-pos].endswith('AG'):
            possible.append(i)
    return possible


def find_all_splice_AC(seq, offset=6):
    # check reference end for splice acceptor, might be off by a few bp when small exons missed
    possible = []
    for i in range(0, offset):
        pos = offset-i
        if seq[:-pos].endswith('AC'):
            possible.append(i)
    return possible


def find_all_splice_CT(seq, offset=6):
    # check reference end for splice acceptor, might be off by a few bp when small exons missed
    possible = []
    for i in range(0, offset):
        pos = offset-i
        if seq[pos:pos+2].startswith('CT'):
            possible.append(i)
    return possible


def filter4_splice_GT(seq, align):
    # given sequence and edlib object, loop through and filter for GT splice sites
    filter = []
    for x in align['locations']:
        if seq[x[0]:x[1]+3].endswith('GT'):
            filter.append(x)
    align['locations'] = filter
    return align


def filter4_splice_AC(seq, align):
    # given sequence and edlib object, loop through and filter for AC splice sites
    filter = []
    for x in align['locations']:
        if seq[x[0]-2:x[1]+1].startswith('AC'):
            filter.append(x)
    align['locations'] = filter
    return align


def filter4_splice_AG(seq, align):
    filter = []
    for x in align['locations']:
        if seq[x[0]-2:x[1]+1].startswith('AG'):
            filter.append(x)
    align['locations'] = filter
    return align


def filter4_splice_CT(seq, align):
    filter = []
    for x in align['locations']:
        if seq[x[0]:x[1]+1+2].endswith('CT'):
            filter.append(x)
    align['locations'] = filter
    return align


def left_update_paf_plus(paf, align, slide, offset):
    if len(align['locations']) > 0:
        old_st = paf[7]
        lp = align['locations'][-1]
        new_st = lp[0]+offset
        paf[2] = 0
        paf[7] = new_st
        old_cs = paf[-1].replace('cs:Z:', '')
        old_cs_tup = cs2tuples(old_cs)
        if slide > 0:
            old_value = int(old_cs_tup[0][1])
            new_value = old_value-slide
            if new_value < 1:
                if old_cs_tup[1][0] == '*':
                    if old_value + (len(old_cs_tup[1][1])/2) == slide:
                        del old_cs_tup[1]
                elif old_cs_tup[1][0] == ':':
                    old_value2 = int(old_cs_tup[1][1])
                    new_value2 = old_value2 + new_value
                    old_cs_tup[1] = (':', new_value2)
                del old_cs_tup[0]
            else:
                old_cs_tup[0] = (':', new_value)
        new_cs_str = ''
        for x in old_cs_tup:
            new_cs_str += '{}{}'.format(x[0], x[1])
        intron_len = (old_st+slide) - (lp[1]+1+offset)
        assert intron_len > 0, "Introns cannot be less than zero"
        new_cs = 'cs:Z::{}~gt{}ag{}'.format(align['cigar'].strip('='), intron_len, new_cs_str)
        try:
            matches, bases = cs2matches(new_cs)
        except ValueError:
            print(paf)
            sys.exit(1)
        paf[9] = matches
        paf[10] = bases
        paf[-1] = new_cs
        return paf
    else:
        return paf


def left_update_paf_minus(paf, align, slide, offset):
    if len(align['locations']) > 0:
        old_en = paf[8]
        lp = align['locations'][0]
        new_en = lp[1]+1+offset
        paf[2] = 0
        paf[8] = new_en
        old_cs = paf[-1].replace('cs:Z:', '')
        old_cs_tup = cs2tuples(old_cs)
        if slide > 0:
            old_value = int(old_cs_tup[-1][1])
            new_value = old_value-slide
            if new_value < 1:
                if old_cs_tup[-2][0] == '*':
                    if old_value + (len(old_cs_tup[-2][1])/2) == slide:
                        del old_cs_tup[-2]
                elif old_cs_tup[-2][0] == ':':
                    old_value2 = int(old_cs_tup[-2][1])
                    new_value2 = old_value2 + new_value
                    old_cs_tup[-2] = (':', new_value2)
                del old_cs_tup[-1]
            else:
                old_cs_tup[-1] = (':', new_value)
        new_cs_str = ''
        for x in old_cs_tup:
            new_cs_str += '{}{}'.format(x[0], x[1])
        intron_len = (lp[0]+offset) - (old_en-slide)
        assert intron_len > 0, "Introns cannot be less than zero"
        new_cs = 'cs:Z:{}~ct{}ac:{}'.format(new_cs_str, intron_len, align['cigar'].strip('='))
        try:
            matches, bases = cs2matches(new_cs)
        except ValueError:
            print(paf)
            sys.exit(1)
        paf[9] = matches
        paf[10] = bases
        paf[-1] = new_cs
        return paf
    else:
        return paf


def right_update_paf_plus(paf, align, slide, offset):
    if len(align['locations']) > 0:
        old_en = paf[8]
        lp = align['locations'][0]
        new_en = lp[1]+1+offset
        paf[3] = paf[1]
        paf[8] = new_en
        old_cs = paf[-1].replace('cs:Z:', '')
        old_cs_tup = cs2tuples(old_cs)
        if slide > 0:
            old_value = int(old_cs_tup[-1][1])
            new_value = old_value-slide
            if new_value < 1:  # this could happen on bad alignment, but can only go a few bp, check next value
                if old_cs_tup[-2][0] == '*':
                    if old_value + (len(old_cs_tup[-2][1])/2) == slide:
                        del old_cs_tup[-2]
                elif old_cs_tup[-2][0] == ':':
                    old_value2 = int(old_cs_tup[-2][1])
                    new_value2 = old_value2 + new_value
                    old_cs_tup[-2] = (':', new_value2)
                del old_cs_tup[-1]
            else:
                old_cs_tup[-1] = (':', new_value)
        new_cs_str = ''
        for x in old_cs_tup:
            new_cs_str += '{}{}'.format(x[0], x[1])
        intron_len = (lp[0]+offset) - (old_en-slide)
        assert intron_len > 0, "Introns cannot be less than zero"
        new_cs = 'cs:Z:{}~gt{}ag:{}'.format(old_cs, intron_len, align['cigar'].strip('='))
        try:
            matches, bases = cs2matches(new_cs)
        except ValueError:
            print(paf)
            sys.exit(1)
        paf[9] = matches
        paf[10] = bases
        paf[-1] = new_cs
        return paf
    else:
        return paf


def right_update_paf_minus(paf, align, slide, offset):
    if len(align['locations']) > 0:
        old_st = paf[7]
        lp = align['locations'][-1]
        new_st = lp[0]+offset
        paf[3] = paf[1]
        paf[7] = new_st
        old_cs = paf[-1].replace('cs:Z:', '')
        old_cs_tup = cs2tuples(old_cs)
        if slide > 0:
            old_value = int(old_cs_tup[0][1])
            new_value = old_value-slide
            if new_value < 1:
                if old_cs_tup[1][0] == '*':
                    if old_value + (len(old_cs_tup[1][1])/2) == slide:
                        del old_cs_tup[1]
                elif old_cs_tup[1][0] == ':':
                    old_value2 = int(old_cs_tup[1][1])
                    new_value2 = old_value2 + new_value
                    old_cs_tup[1] = (':', new_value2)
                del old_cs_tup[0]
            else:
                old_cs_tup[0] = (':', new_value)
        new_cs_str = ''
        for x in old_cs_tup:
            new_cs_str += '{}{}'.format(x[0], x[1])
        intron_len = (old_st-slide) - (lp[1]+1+offset)
        assert intron_len > 0, "Introns cannot be less than zero"
        new_cs = 'cs:Z::{}~ct{}ac{}'.format(align['cigar'].strip('='), intron_len, new_cs_str)
        try:
            matches, bases = cs2matches(new_cs)
        except ValueError:
            print(paf)
            print(new_cs)
            sys.exit(1)
        paf[9] = matches
        paf[10] = bases
        paf[-1] = new_cs
        return paf
    else:
        return paf


def cs2matches(cs):
    # return num matching and number of bases
    tups = cs2tuples(cs.replace('cs:Z:', ''))
    matches = 0
    bases = 0
    for t in tups:
        if t[0] == ':':
            matches += int(t[1])
            bases += int(t[1])
        elif t[0] in ['+', '-']:
            bases += len(t[1])
    return matches, bases


def cs2tuples(cs, separators=[':', '*', '+', '-', '~']):
    # separators is an array of strings that are being used to split the the string.
    tmpList = []
    i = 0
    while i < len(cs):
        theSeparator = ""
        for current in separators:
            if current == cs[i:i+len(current)]:
                theSeparator = current
        if theSeparator != "":
            tmpList += [theSeparator]
            i = i + len(theSeparator)
        else:
            if tmpList == []:
                tmpList = [""]
            if(tmpList[-1] in separators):
                tmpList += [""]
            tmpList[-1] += cs[i]
            i += 1
    tupList = list(zip(tmpList[::2], tmpList[1::2]))
    return tupList


def left_plus_best_align(paf, refseq, contig, refstart, queryseq, querystart, offset=6, maxlen=500):
    keep =[]
    if refstart-maxlen < 0:
        r_st = 0
    else:
        r_st = refstart-maxlen
    ref = refseq.seq(contig, r_st, refstart+offset)
    slides = find_all_splice_AG(ref, offset=offset)
    for slide in slides:  # look through possible splice sites
        query = queryseq[0:querystart+slide]
        align = edlib.align(query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)
        align = filter4_splice_GT(ref, align)
        if len(align['locations']) < 1:
            continue
        else:
            keep.append((align, slide))
    if len(keep) > 0:
        a, s = keep[0]
        paf = left_update_paf_plus(paf, a, s, refstart-maxlen)
    return paf


def left_minus_best_align(paf, refseq, contig, refend, queryseq, querystart, offset=6, maxlen=500):
    keep =[]
    if refend+maxlen > len(refseq.seq(contig)):
        r_en = len(refseq.seq(contig))
    else:
        r_en = refend+maxlen
    ref = refseq.seq(contig, refend-offset, r_en)
    slides = find_all_splice_CT(ref, offset=offset)
    for slide in slides:  # look through possible splice sites
        query = mp.revcomp(queryseq[0:querystart+slide])
        align = edlib.align(query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)
        align = filter4_splice_AC(ref, align)
        if len(align['locations']) < 1:
            continue
        else:
            keep.append((align, slide))
    if len(keep) > 0:
        a, s = keep[0]
        paf = left_update_paf_minus(paf, a, s, refend-offset)
    return paf


def right_plus_best_align(paf, refseq, contig, refend, queryseq, queryend, offset=6, maxlen=500):
    keep =[]
    if refend+maxlen > len(refseq.seq(contig)):
        r_en = len(refseq.seq(contig))
    else:
        r_en = refend+maxlen
    ref = refseq.seq(contig, refend-offset, r_en)
    slides = find_all_splice_GT(ref, offset=offset)
    for slide in slides:  # look through possible splice sites
        query = queryseq[queryend-slide:]
        align = edlib.align(query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)
        align = filter4_splice_AG(ref, align)
        if len(align['locations']) < 1:
            continue
        else:
            keep.append((align, slide))
    if len(keep) > 0:
        a, s = keep[0]
        paf = right_update_paf_plus(paf, a, s, refend-offset)
    return paf


def right_minus_best_align(paf, refseq, contig, refstart, queryseq, queryend, offset=6, maxlen=500):
    keep = []
    if refstart-maxlen < 0:
        r_st = 0
    else:
        r_st = refstart-maxlen
    ref = refseq.seq(contig, r_st, refstart+offset)
    slides = find_all_splice_AC(ref, offset=offset)
    for slide in slides:  # look through possible splice sites
        query = mp.revcomp(queryseq[queryend-slide:])
        align = edlib.align(query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)
        align = filter4_splice_CT(ref, align)
        if len(align['locations']) < 1:
            continue
        else:
            keep.append((align, slide))
    if len(keep) > 0:
        a, s = keep[0]
        paf = right_update_paf_minus(paf, a, s, refstart-maxlen)
    return paf


def splice_aligner(reference, query, threads=3, min_mapq=1, max_intron=500):
    """
    # splicing plus strand, GT-AG junctions
    ATG GGC | GTX  ------- XAG | GCC TAA

    # splicing rev strand, CT-AC junctions
    TTA GGC | CTX ------- XAC GCC CAT

    """
    results = []
    stats = {'n': 0, 'low-mapq': 0, 'refine-left': 0, 'refine-right': 0}
    # run minimap2 splice alignment and refine alignment with edlib
    # load reference into index for search and index access
    ref_idx = mp.Aligner(reference, preset='splice', n_threads=threads)
    if not ref_idx: raise Exception("ERROR: failed to load/build index of {}".format(reference))
    # now go through query aligning to ref
    for name, seq, qual in mp.fastx_read(query):
        for h in ref_idx.map(seq, cs=True):
            stats['n'] += 1
            if h.mapq < min_mapq:
                stats['low-mapq'] += 1
                continue
            # construct data for paf format
            if h.strand > 0: strand = '+'
            elif h.strand < 0: strand = '-'
            else: strand = '?'
            nm = 'NM:i:{}'.format(h.NM)
            if h.is_primary != 0: tp = 'tp:A:P'
            else: tp = 'tp:A:S'
            if h.trans_strand > 0: ts = 'ts:A:+'
            elif h.trans_strand < 0: ts = 'ts:A:-'
            else: ts = 'ts:A:.'
            cs = 'cs:Z:{}'.format(h.cs)
            # store paf format as a list
            paf = [name, len(seq), h.q_st, h.q_en, strand, h.ctg,
                   h.ctg_len, h.r_st, h.r_en, h.mlen, h.blen, h.mapq,
                   tp, ts, nm, cs]
            if h.q_st > 0: # refine the alignment for first exon or left side of query
                if strand == '+':
                    paf = left_plus_best_align(paf, ref_idx, h.ctg, h.r_st, seq, h.q_st, offset=6, maxlen=max_intron)
                else:  # negative/crick strand
                    paf = left_minus_best_align(paf, ref_idx, h.ctg, h.r_en, seq, h.q_st, offset=6, maxlen=max_intron)
                if paf[-1] != cs:
                    stats['refine-left'] += 1
            if len(seq) > h.q_en:  # refine alignment for last exon or right side of query
                if strand == '+':
                    paf = right_plus_best_align(paf, ref_idx, h.ctg, h.r_en, seq, h.q_en, offset=6, maxlen=max_intron)
                else:
                    paf = right_minus_best_align(paf, ref_idx, h.ctg, h.r_st, seq, h.q_en, offset=6, maxlen=max_intron)
                if paf[-1] != cs:
                    stats['refine-right'] += 1
            results.append(paf)
    return results, stats


def aligner(reference, query, output=False, threads=3, min_mapq=1, max_intron=500, debug=False):
    # wrapper for spliced alignment for script to write to file
    if output:
        outfile = zopen(output, mode='w')
    else:
        outfile = sys.stdout
    results, stats = splice_aligner(reference, query, threads=threads, min_mapq=min_mapq, max_intron=max_intron)
    if debug:
        sys.stderr.write("Generated {} alignments: {} required 5' refinement, {} required 3' refinement, {} dropped low score\n".format(
            stats['n'], stats['refine-left'], stats['refine-right'], stats['low-mapq']))
    for r in results:
        outfile.write('{}\n'.format('\t'.join([str(x) for x in r])))


def align(args):
    check_inputs([args.reference, args.query])
    aligner(args.reference, args.query, output=args.out, threads=args.threads,
            min_mapq=args.min_mapq, max_intron=args.max_intron, debug=args.debug)





