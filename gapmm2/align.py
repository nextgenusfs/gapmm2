import sys

import edlib
import mappy as mp
from natsort import natsorted

from .utils import check_inputs, zopen

degenNuc = [
    ("R", "A"),
    ("R", "G"),
    ("M", "A"),
    ("M", "C"),
    ("W", "A"),
    ("W", "T"),
    ("S", "C"),
    ("S", "G"),
    ("Y", "C"),
    ("Y", "T"),
    ("K", "G"),
    ("K", "T"),
    ("V", "A"),
    ("V", "C"),
    ("V", "G"),
    ("H", "A"),
    ("H", "C"),
    ("H", "T"),
    ("D", "A"),
    ("D", "G"),
    ("D", "T"),
    ("B", "C"),
    ("B", "G"),
    ("B", "T"),
    ("N", "G"),
    ("N", "A"),
    ("N", "T"),
    ("N", "C"),
    ("X", "G"),
    ("X", "A"),
    ("X", "T"),
    ("X", "C"),
]


def find_all_splice_GT(seq, offset=6):
    """
    Find all potential GT splice donor sites in a sequence.

    This function checks the reference end for splice donor sites, which might be off by a few
    base pairs when small exons are missed.

    Args:
        seq (str): The DNA sequence to search.
        offset (int, optional): The maximum offset to search. Defaults to 6.

    Returns:
        list: A list of positions where GT splice donor sites were found.
    """
    possible = []
    for i in range(0, offset):
        pos = offset - i
        if seq[pos : pos + 2].startswith("GT"):
            possible.append(i)
    return possible


def find_all_splice_AG(seq, offset=6):
    """
    Find all potential AG splice acceptor sites in a sequence.

    This function checks the reference end for splice acceptor sites, which might be off by a few
    base pairs when small exons are missed.

    Args:
        seq (str): The DNA sequence to search.
        offset (int, optional): The maximum offset to search. Defaults to 6.

    Returns:
        list: A list of positions where AG splice acceptor sites were found.
    """
    possible = []
    for i in range(0, offset):
        pos = offset - i
        if seq[:-pos].endswith("AG"):
            possible.append(i)
    return possible


def find_all_splice_AC(seq, offset=6):
    """
    Find all potential AC splice acceptor sites in a sequence.

    This function checks the reference end for splice acceptor sites, which might be off by a few
    base pairs when small exons are missed.

    Args:
        seq (str): The DNA sequence to search.
        offset (int, optional): The maximum offset to search. Defaults to 6.

    Returns:
        list: A list of positions where AC splice acceptor sites were found.
    """
    possible = []
    for i in range(0, offset):
        pos = offset - i
        if seq[:-pos].endswith("AC"):
            possible.append(i)
    return possible


def find_all_splice_CT(seq, offset=6):
    """
    Find all potential CT splice donor sites in a sequence.

    This function checks the reference end for splice donor sites, which might be off by a few
    base pairs when small exons are missed.

    Args:
        seq (str): The DNA sequence to search.
        offset (int, optional): The maximum offset to search. Defaults to 6.

    Returns:
        list: A list of positions where CT splice donor sites were found.
    """
    possible = []
    for i in range(0, offset):
        pos = offset - i
        if seq[pos : pos + 2].startswith("CT"):
            possible.append(i)
    return possible


def filter4_splice_GT(seq, align):
    """
    Filter alignment locations for GT splice donor sites.

    Given a sequence and an edlib alignment object, this function filters the alignment
    locations to keep only those that end with a GT splice donor site.

    Args:
        seq (str): The DNA sequence to search.
        align (dict): An edlib alignment object with a 'locations' key.

    Returns:
        dict: The filtered edlib alignment object.
    """
    filter = []
    for x in align["locations"]:
        if seq[x[0] : x[1] + 3].endswith("GT"):
            filter.append(x)
    align["locations"] = filter
    return align


def filter4_splice_AC(seq, align):
    """
    Filter alignment locations for AC splice acceptor sites.

    Given a sequence and an edlib alignment object, this function filters the alignment
    locations to keep only those that start with an AC splice acceptor site.

    Args:
        seq (str): The DNA sequence to search.
        align (dict): An edlib alignment object with a 'locations' key.

    Returns:
        dict: The filtered edlib alignment object.
    """
    filter = []
    for x in align["locations"]:
        if seq[x[0] - 2 : x[1] + 1].startswith("AC"):
            filter.append(x)
    align["locations"] = filter
    return align


def filter4_splice_AG(seq, align):
    """
    Filter alignment locations for AG splice acceptor sites.

    Given a sequence and an edlib alignment object, this function filters the alignment
    locations to keep only those that start with an AG splice acceptor site.

    Args:
        seq (str): The DNA sequence to search.
        align (dict): An edlib alignment object with a 'locations' key.

    Returns:
        dict: The filtered edlib alignment object.
    """
    filter = []
    for x in align["locations"]:
        if seq[x[0] - 2 : x[1] + 1].startswith("AG"):
            filter.append(x)
    align["locations"] = filter
    return align


def filter4_splice_CT(seq, align):
    """
    Filter alignment locations for CT splice donor sites.

    Given a sequence and an edlib alignment object, this function filters the alignment
    locations to keep only those that end with a CT splice donor site.

    Args:
        seq (str): The DNA sequence to search.
        align (dict): An edlib alignment object with a 'locations' key.

    Returns:
        dict: The filtered edlib alignment object.
    """
    filter = []
    for x in align["locations"]:
        if seq[x[0] : x[1] + 1 + 2].endswith("CT"):
            filter.append(x)
    align["locations"] = filter
    return align


def left_update_paf_plus(paf, align, slide, offset):
    if len(align["locations"]) > 0:
        old_st = paf[7]
        lp = align["locations"][-1]
        new_st = lp[0] + offset
        paf[2] = 0
        paf[7] = new_st
        old_cs = paf[-1].replace("cs:Z:", "")
        old_cs_tup = cs2tuples(old_cs)
        if slide > 0:
            old_value = int(old_cs_tup[0][1])
            new_value = old_value - slide
            if new_value < 1:
                if old_cs_tup[1][0] == "*":
                    if old_value + (len(old_cs_tup[1][1]) / 2) == slide:
                        del old_cs_tup[1]
                elif old_cs_tup[1][0] == ":":
                    old_value2 = int(old_cs_tup[1][1])
                    new_value2 = old_value2 + new_value
                    old_cs_tup[1] = (":", new_value2)
                del old_cs_tup[0]
            else:
                old_cs_tup[0] = (":", new_value)
        new_cs_str = ""
        for x in old_cs_tup:
            new_cs_str += "{}{}".format(x[0], x[1])
        intron_len = (old_st + slide) - (lp[1] + 1 + offset)
        try:
            assert intron_len > 0, "Introns cannot be less than zero"
        except AssertionError:
            return paf
        new_cs = "cs:Z::{}~gt{}ag{}".format(align["cigar"].strip("="), intron_len, new_cs_str)
        try:
            matches, bases = cs2matches(new_cs)
        except ValueError:
            print(paf)
            raise SystemExit(1)
        paf[9] = matches
        paf[10] = bases
        paf[-1] = new_cs
        return paf
    else:
        return paf


def left_update_paf_minus(paf, align, slide, offset):
    if len(align["locations"]) > 0:
        old_en = paf[8]
        lp = align["locations"][0]
        new_en = lp[1] + 1 + offset
        paf[2] = 0
        paf[8] = new_en
        old_cs = paf[-1].replace("cs:Z:", "")
        old_cs_tup = cs2tuples(old_cs)
        if slide > 0:
            old_value = int(old_cs_tup[-1][1])
            new_value = old_value - slide
            if new_value < 1:
                if old_cs_tup[-2][0] == "*":
                    if old_value + (len(old_cs_tup[-2][1]) / 2) == slide:
                        del old_cs_tup[-2]
                elif old_cs_tup[-2][0] == ":":
                    old_value2 = int(old_cs_tup[-2][1])
                    new_value2 = old_value2 + new_value
                    old_cs_tup[-2] = (":", new_value2)
                del old_cs_tup[-1]
            else:
                old_cs_tup[-1] = (":", new_value)
        new_cs_str = ""
        for x in old_cs_tup:
            new_cs_str += "{}{}".format(x[0], x[1])
        intron_len = (lp[0] + offset) - (old_en - slide)
        try:
            assert intron_len > 0, "Introns cannot be less than zero"
        except AssertionError:
            return paf
        new_cs = "cs:Z:{}~ct{}ac:{}".format(new_cs_str, intron_len, align["cigar"].strip("="))
        try:
            matches, bases = cs2matches(new_cs)
        except ValueError:
            print(paf)
            raise SystemExit(1)
        paf[9] = matches
        paf[10] = bases
        paf[-1] = new_cs
        return paf
    else:
        return paf


def right_update_paf_plus(paf, align, slide, offset):
    if len(align["locations"]) > 0:
        old_en = paf[8]
        lp = align["locations"][0]
        new_en = lp[1] + 1 + offset
        paf[3] = paf[1]
        paf[8] = new_en
        old_cs = paf[-1].replace("cs:Z:", "")
        old_cs_tup = cs2tuples(old_cs)
        if slide > 0:
            old_value = int(old_cs_tup[-1][1])
            new_value = old_value - slide
            if (
                new_value < 1
            ):  # this could happen on bad alignment, but can only go a few bp, check next value
                if old_cs_tup[-2][0] == "*":
                    if old_value + (len(old_cs_tup[-2][1]) / 2) == slide:
                        del old_cs_tup[-2]
                elif old_cs_tup[-2][0] == ":":
                    old_value2 = int(old_cs_tup[-2][1])
                    new_value2 = old_value2 + new_value
                    old_cs_tup[-2] = (":", new_value2)
                del old_cs_tup[-1]
            else:
                old_cs_tup[-1] = (":", new_value)
        new_cs_str = ""
        for x in old_cs_tup:
            new_cs_str += "{}{}".format(x[0], x[1])
        intron_len = (lp[0] + offset) - (old_en - slide)
        try:
            assert intron_len > 0, "Introns cannot be less than zero"
        except AssertionError:
            return paf
        new_cs = "cs:Z:{}~gt{}ag:{}".format(old_cs, intron_len, align["cigar"].strip("="))
        try:
            matches, bases = cs2matches(new_cs)
        except ValueError:
            print(paf)
            raise SystemExit(1)
        paf[9] = matches
        paf[10] = bases
        paf[-1] = new_cs
        return paf
    else:
        return paf


def right_update_paf_minus(paf, align, slide, offset):
    if len(align["locations"]) > 0:
        old_st = paf[7]
        lp = align["locations"][-1]
        new_st = lp[0] + offset
        paf[3] = paf[1]
        paf[7] = new_st
        old_cs = paf[-1].replace("cs:Z:", "")
        old_cs_tup = cs2tuples(old_cs)
        if slide > 0:
            old_value = int(old_cs_tup[0][1])
            new_value = old_value - slide
            if new_value < 1:
                if old_cs_tup[1][0] == "*":
                    if old_value + (len(old_cs_tup[1][1]) / 2) == slide:
                        del old_cs_tup[1]
                elif old_cs_tup[1][0] == ":":
                    old_value2 = int(old_cs_tup[1][1])
                    new_value2 = old_value2 + new_value
                    old_cs_tup[1] = (":", new_value2)
                del old_cs_tup[0]
            else:
                old_cs_tup[0] = (":", new_value)
        new_cs_str = ""
        for x in old_cs_tup:
            new_cs_str += "{}{}".format(x[0], x[1])
        intron_len = (old_st - slide) - (lp[1] + 1 + offset)
        try:
            assert intron_len > 0, "Introns cannot be less than zero"
        except AssertionError:
            return paf
        new_cs = "cs:Z::{}~ct{}ac{}".format(align["cigar"].strip("="), intron_len, new_cs_str)
        try:
            matches, bases = cs2matches(new_cs)
        except ValueError:
            print(paf)
            print(new_cs)
            raise SystemExit(1)
        paf[9] = matches
        paf[10] = bases
        paf[-1] = new_cs
        return paf
    else:
        return paf


def cs2matches(cs):
    # return num matching and number of bases
    tups = cs2tuples(cs.replace("cs:Z:", ""))
    matches = 0
    bases = 0
    for t in tups:
        if t[0] == ":":
            matches += int(t[1])
            bases += int(t[1])
        elif t[0] in ["+", "-"]:
            bases += len(t[1])
    return matches, bases


def cs2tuples(cs, separators=[":", "*", "+", "-", "~"]):
    # separators is an array of strings that are being used to split the the string.
    tmpList = []
    i = 0
    while i < len(cs):
        theSeparator = ""
        for current in separators:
            if current == cs[i : i + len(current)]:
                theSeparator = current
        if theSeparator != "":
            tmpList += [theSeparator]
            i = i + len(theSeparator)
        else:
            if tmpList == []:
                tmpList = [""]
            if tmpList[-1] in separators:
                tmpList += [""]
            tmpList[-1] += cs[i]
            i += 1
    tupList = list(zip(tmpList[::2], tmpList[1::2]))
    return tupList


def left_plus_best_align(paf, refseq, contig, refstart, queryseq, querystart, offset=6, maxlen=500):
    keep = []
    if refstart - maxlen < 0:
        r_st = 0
    else:
        r_st = refstart - maxlen
    ref = refseq.seq(contig, r_st, refstart + offset)
    slides = find_all_splice_AG(ref, offset=offset)
    for slide in slides:  # look through possible splice sites
        query = queryseq[0 : querystart + slide]
        align = edlib.align(query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)
        align = filter4_splice_GT(ref, align)
        if len(align["locations"]) < 1:
            continue
        else:
            keep.append((align, slide))
    if len(keep) > 0:
        a, s = keep[0]
        paf = left_update_paf_plus(paf, a, s, refstart - maxlen)
    return paf


def left_minus_best_align(paf, refseq, contig, refend, queryseq, querystart, offset=6, maxlen=500):
    keep = []
    if refend + maxlen > len(refseq.seq(contig)):
        r_en = len(refseq.seq(contig))
    else:
        r_en = refend + maxlen
    ref = refseq.seq(contig, refend - offset, r_en)
    slides = find_all_splice_CT(ref, offset=offset)
    for slide in slides:  # look through possible splice sites
        query = mp.revcomp(queryseq[0 : querystart + slide])
        align = edlib.align(query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)
        align = filter4_splice_AC(ref, align)
        if len(align["locations"]) < 1:
            continue
        else:
            keep.append((align, slide))
    if len(keep) > 0:
        a, s = keep[0]
        paf = left_update_paf_minus(paf, a, s, refend - offset)
    return paf


def right_plus_best_align(paf, refseq, contig, refend, queryseq, queryend, offset=6, maxlen=500):
    keep = []
    if refend + maxlen > len(refseq.seq(contig)):
        r_en = len(refseq.seq(contig))
    else:
        r_en = refend + maxlen
    ref = refseq.seq(contig, refend - offset, r_en)
    slides = find_all_splice_GT(ref, offset=offset)
    for slide in slides:  # look through possible splice sites
        query = queryseq[queryend - slide :]
        align = edlib.align(query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)
        align = filter4_splice_AG(ref, align)
        if len(align["locations"]) < 1:
            continue
        else:
            keep.append((align, slide))
    if len(keep) > 0:
        a, s = keep[0]
        paf = right_update_paf_plus(paf, a, s, refend - offset)
    return paf


def right_minus_best_align(paf, refseq, contig, refstart, queryseq, queryend, offset=6, maxlen=500):
    keep = []
    if refstart - maxlen < 0:
        r_st = 0
    else:
        r_st = refstart - maxlen
    ref = refseq.seq(contig, r_st, refstart + offset)
    slides = find_all_splice_AC(ref, offset=offset)
    for slide in slides:  # look through possible splice sites
        query = mp.revcomp(queryseq[queryend - slide :])
        align = edlib.align(query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)
        align = filter4_splice_CT(ref, align)
        if len(align["locations"]) < 1:
            continue
        else:
            keep.append((align, slide))
    if len(keep) > 0:
        a, s = keep[0]
        paf = right_update_paf_minus(paf, a, s, refstart - maxlen)
    return paf


def splice_aligner(reference, query, threads=3, min_mapq=1, max_intron=500):
    """
    Perform spliced alignment of transcripts to a genome using minimap2 and refine terminal alignments with edlib.

    This function aligns transcripts to a genome using minimap2 with splice options, and then refines
    the terminal alignments using edlib to improve the alignment of terminal exons that might be missed
    by minimap2. The function handles both forward and reverse strand alignments, recognizing the
    appropriate splice junctions for each strand.

    Splice junctions:
    - Plus strand: GT-AG junctions (ATG GGC | GTX  ------- XAG | GCC TAA)
    - Reverse strand: CT-AC junctions (TTA GGC | CTX ------- XAC GCC CAT)

    Args:
        reference (str): Path to the reference genome FASTA file.
        query (str): Path to the query transcripts FASTA or FASTQ file.
        threads (int, optional): Number of threads to use with minimap2. Defaults to 3.
        min_mapq (int, optional): Minimum mapping quality value to keep an alignment. Defaults to 1.
        max_intron (int, optional): Maximum intron length, controls terminal search space. Defaults to 500.

    Returns:
        tuple: A tuple containing:
            - list: A list of lists, where each inner list contains PAF formatted data for an alignment.
            - dict: A dictionary with statistics about the alignments, including:
                - 'n': Total number of alignments.
                - 'low-mapq': Number of alignments dropped due to low mapping quality.
                - 'refine-left': Number of alignments where the left side was refined.
                - 'refine-right': Number of alignments where the right side was refined.
    """
    results = []
    stats = {"n": 0, "low-mapq": 0, "refine-left": 0, "refine-right": 0}
    # run minimap2 splice alignment and refine alignment with edlib
    # load reference into index for search and index access
    ref_idx = mp.Aligner(reference, preset="splice", n_threads=threads)
    if not ref_idx:
        raise Exception("ERROR: failed to load/build index of {}".format(reference))
    # now go through query aligning to ref
    for name, seq, _ in mp.fastx_read(query):
        for h in ref_idx.map(seq, cs=True):
            stats["n"] += 1
            if h.mapq < min_mapq:
                stats["low-mapq"] += 1
                continue
            # construct data for paf format
            if h.strand > 0:
                strand = "+"
            elif h.strand < 0:
                strand = "-"
            else:
                strand = "?"
            nm = "NM:i:{}".format(h.NM)
            if h.is_primary != 0:
                tp = "tp:A:P"
            else:
                tp = "tp:A:S"
            if h.trans_strand > 0:
                ts = "ts:A:+"
            elif h.trans_strand < 0:
                ts = "ts:A:-"
            else:
                ts = "ts:A:."
            cs = "cs:Z:{}".format(h.cs)
            # store paf format as a list
            paf = [
                name,
                len(seq),
                h.q_st,
                h.q_en,
                strand,
                h.ctg,
                h.ctg_len,
                h.r_st,
                h.r_en,
                h.mlen,
                h.blen,
                h.mapq,
                tp,
                ts,
                nm,
                cs,
            ]
            if h.q_st > 0:  # refine the alignment for first exon or left side of query
                if strand == "+":
                    paf = left_plus_best_align(
                        paf,
                        ref_idx,
                        h.ctg,
                        h.r_st,
                        seq,
                        h.q_st,
                        offset=6,
                        maxlen=max_intron,
                    )
                else:  # negative/crick strand
                    paf = left_minus_best_align(
                        paf,
                        ref_idx,
                        h.ctg,
                        h.r_en,
                        seq,
                        h.q_st,
                        offset=6,
                        maxlen=max_intron,
                    )
                if paf[-1] != cs:
                    stats["refine-left"] += 1
            if len(seq) > h.q_en:  # refine alignment for last exon or right side of query
                if strand == "+":
                    paf = right_plus_best_align(
                        paf,
                        ref_idx,
                        h.ctg,
                        h.r_en,
                        seq,
                        h.q_en,
                        offset=6,
                        maxlen=max_intron,
                    )
                else:
                    paf = right_minus_best_align(
                        paf,
                        ref_idx,
                        h.ctg,
                        h.r_st,
                        seq,
                        h.q_en,
                        offset=6,
                        maxlen=max_intron,
                    )
                if paf[-1] != cs:
                    stats["refine-right"] += 1
            results.append(paf)
    return results, stats


def cs2coords(
    start,
    qstart,
    length,
    strand,
    cs,
    offset=1,
    splice_donor=["gt", "at"],
    splice_acceptor=["ag", "ac"],
):
    """
    Convert a CIGAR string (cs) from minimap2 to genomic coordinates.

    This function parses the CIGAR string from minimap2 and converts it to genomic coordinates,
    identifying exons, introns, and other alignment features. It handles both forward and reverse
    strand alignments, and can recognize proper splice sites.

    From minimap2 manual, the cs flag definitions are:
    - Op: =, Regex: [ACGTN]+, Description: Identical sequence (long form)
    - Op: :, Regex: [0-9]+, Description: Identical sequence length
    - Op: *, Regex: [acgtn][acgtn], Description: Substitution: ref to query
    - Op: +, Regex: [acgtn]+, Description: Insertion to the reference
    - Op: -, Regex: [acgtn]+, Description: Deletion from the reference
    - Op: ~, Regex: [acgtn]{2}[0-9]+[acgtn]{2}, Description: Intron length and splice signal

    Args:
        start (int): The start position in the reference genome.
        qstart (int): The start position in the query sequence.
        length (int): The length of the query sequence.
        strand (str): The strand of the alignment ("+" or "-").
        cs (str): The CIGAR string from minimap2.
        offset (int, optional): Offset to add to the exon coordinates. Defaults to 1.
        splice_donor (list, optional): List of splice donor sequences. Defaults to ["gt", "at"].
        splice_acceptor (list, optional): List of splice acceptor sequences. Defaults to ["ag", "ac"].

    Returns:
        tuple: A tuple containing:
            - list: A list of exon coordinates as tuples (start, end).
            - list: A list of query coordinates as tuples (start, end).
            - int: Number of mismatches in the alignment.
            - int: Number of gaps in the alignment.
            - bool: Whether the alignment has proper splice sites.
    """
    cs = cs.replace("cs:Z:", "")
    ProperSplice = True
    exons = [int(start)]
    position = int(start)
    query = [int(qstart)]
    querypos = 0
    num_exons = 1
    gaps = 0
    mismatches = 0
    indels = []
    if strand == "+":
        sp_donor = splice_donor
        sp_acceptor = splice_acceptor
        sort_orientation = False
    elif strand == "-":
        # rev comp and swap donor/acceptor
        sp_donor = [mp.revcomp(x).lower() for x in splice_acceptor]
        sp_acceptor = [mp.revcomp(x).lower() for x in splice_donor]
        sort_orientation = True
    for s, value in cs2tuples(cs):
        if s == ":":
            position += int(value)
            querypos += int(value)
            indels.append(0)
        elif s == "-":
            gaps += 1
            position += len(value)
            querypos += len(value)
            indels.append(-len(value))
        elif s == "+":
            gaps += 1
            position += len(value)
            querypos += len(value)
            indels.append(len(value))
        elif s == "~":
            if value.startswith(tuple(sp_donor)) and value.endswith(tuple(sp_acceptor)):
                ProperSplice = True
            else:
                ProperSplice = False
            num_exons += 1
            exons.append(position + indels[-1])
            query.append(querypos)
            query.append(querypos + 1)
            intronLen = int(value[2:-2])
            position += intronLen
            exons.append(position)
            indels.append(0)
        elif s == "*":
            mismatches += len(value) / 2
    # add last Position
    exons.append(position)
    query.append(int(length))
    # convert exon list into list of exon tuples
    exontmp = list(zip(exons[0::2], exons[1::2]))
    queryList = list(zip(query[0::2], query[1::2]))
    exonList = []
    for x in sorted(exontmp, key=lambda tup: tup[0], reverse=sort_orientation):
        exonList.append((x[0] + offset, x[1]))
    return exonList, queryList, mismatches, gaps, ProperSplice


def paf2gff3(paf, output=False, minpident=0):
    """
    Convert PAF format alignments to GFF3 format.

    This function takes a list of PAF format alignments and converts them to GFF3 format,
    writing the output to a file or stdout. It calculates the percent identity for each
    alignment and filters out alignments below a minimum percent identity threshold.

    Args:
        paf (list): A list of lists, where each inner list contains PAF formatted data for an alignment.
        output (str, optional): Path to the output file. If False, writes to stdout. Defaults to False.
        minpident (float, optional): Minimum percent identity threshold. Alignments below this threshold
            will be filtered out. Defaults to 0.
    """
    if output:
        outfile = zopen(output, mode="w")
    else:
        outfile = sys.stdout
    outfile.write("##gff-version 3\n")
    # paf is a list of lists in paf format
    # paf = [name, len(seq), h.q_st, h.q_en, strand, h.ctg, h.ctg_len, h.r_st, h.r_en, h.mlen, h.blen, h.mapq, tp, ts, nm, cs]
    count = 1
    sorted_paf = natsorted(paf, key=lambda x: (x[5], x[7]))
    for p in sorted_paf:
        exons, queries, mismatches, gaps, _ = cs2coords(p[7], p[2], p[9], p[4], p[-1])
        matches = p[9] - mismatches - gaps
        pident = 100 * (matches / p[9])
        if pident < minpident:
            continue
        for i, exon in enumerate(exons):
            start = exon[0]
            end = exon[1]
            if p[4] == "+":
                qstart = queries[i][0]
                qend = queries[i][1]
            else:
                qstart = p[1] - queries[i][1] + 1
                qend = p[1] - queries[i][0] + 1
            if qstart == 0:
                qstart = 1
            outfile.write(
                "{:}\t{:}\t{:}\t{:}\t{:}\t{:.2f}\t{:}\t{:}\tID=gapmm2_{:};Target={:} {:} {:}\n".format(
                    p[5],
                    "gapmm2",
                    "cDNA_match",
                    start,
                    end,
                    pident,
                    p[4],
                    ".",
                    count,
                    p[0],
                    qstart,
                    qend,
                )
            )
        count += 1


def aligner(
    reference,
    query,
    output=False,
    threads=3,
    out_fmt="paf",
    min_mapq=1,
    max_intron=500,
    debug=False,
):
    """
    Write spliced alignment results to a file or stdout.

    This function is the main entry point for the command-line interface. It performs spliced
    alignment using the splice_aligner function and writes the results to a file or stdout in
    the specified format.

    Args:
        reference (str): Path to the reference genome FASTA file.
        query (str): Path to the query transcripts FASTA or FASTQ file.
        output (str, optional): Path to the output file. If False, writes to stdout. Defaults to False.
        threads (int, optional): Number of threads to use with minimap2. Defaults to 3.
        out_fmt (str, optional): Output format, either "paf" or "gff3". Defaults to "paf".
        min_mapq (int, optional): Minimum mapping quality value to keep an alignment. Defaults to 1.
        max_intron (int, optional): Maximum intron length, controls terminal search space. Defaults to 500.
        debug (bool, optional): Whether to write debug information to stderr. Defaults to False.
    """
    # wrapper for spliced alignment for script to write to file
    if output:
        outfile = zopen(output, mode="w")
    else:
        outfile = sys.stdout
    results, stats = splice_aligner(
        reference, query, threads=threads, min_mapq=min_mapq, max_intron=max_intron
    )
    if debug:
        sys.stderr.write(
            "Generated {} alignments: {} required 5' refinement, {} required 3' refinement, {} dropped low score\n".format(
                stats["n"],
                stats["refine-left"],
                stats["refine-right"],
                stats["low-mapq"],
            )
        )
    if out_fmt == "paf":
        for r in natsorted(results, key=lambda x: (x[5], x[7])):
            outfile.write("{}\n".format("\t".join([str(x) for x in r])))
    elif out_fmt == "gff3":
        paf2gff3(results, output=output)


def align(args):
    """
    Entry point for the command-line interface.

    This function is called from the __main__.py file and is the entry point for the
    command-line interface. It checks the input files and calls the aligner function
    with the appropriate arguments.

    Args:
        args: An argparse.Namespace object containing the command-line arguments.
    """
    check_inputs([args.reference, args.query])
    aligner(
        args.reference,
        args.query,
        output=args.out,
        threads=args.threads,
        min_mapq=args.min_mapq,
        max_intron=args.max_intron,
        out_fmt=args.out_fmt,
        debug=args.debug,
    )
