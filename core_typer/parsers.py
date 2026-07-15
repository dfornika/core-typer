import csv
import logging
import re

logger = logging.getLogger(__name__)

_MD5_HEX_RE = re.compile(r"^[0-9a-f]{32}$")


def parse_template_name(template):
    """
    Split a scheme template name into (locus_id, allele_id, allele_hash).

    Templates are named "<locus_id>_<allele_id>", or, for schemes prepared
    with allele hashing (see scripts/prepare_cgmlst_scheme.py), "<locus_id>_
    <allele_id>_<md5_hash>". allele_hash is None if no hash suffix is
    present. locus_id is everything before the allele_id/hash, so locus
    names containing underscores are still handled correctly.

    :param template: The template name, e.g. "locusA_1" or "locusA_1_<md5>"
    :type template: str
    :return: (locus_id, allele_id, allele_hash)
    :rtype: tuple[str, str, str | None]
    """
    parts = template.split("_")
    if len(parts) >= 3 and _MD5_HEX_RE.match(parts[-1]):
        allele_hash = parts[-1]
        allele_id = parts[-2]
        locus_id = "_".join(parts[:-2])
    else:
        allele_hash = None
        allele_id = parts[-1]
        locus_id = "_".join(parts[:-1])
    return locus_id, allele_id, allele_hash


def parse_kma_result(kma_result_file):
    """
    Parse a kma result file into a dict of lists of dicts.

    :param kma_result_file: The path to the kma result file
    :type kma_result_file: str
    :return: The kma results, indexed by locus_id. Keys of kma results are: locus_id, allele_id, allele_hash, score, template_length, template_identity, template_coverage, query_identity, query_coverage, depth, q_value, p_value
    :rtype: dict[str, list[dict]]
    """
    kma_result_by_locus_id = {}
    int_fields = [
        "score",
        "expected",
        "template_length",
    ]
    float_fields = [
        "template_identity",
        "template_coverage",
        "query_identity",
        "query_coverage",
        "depth",
        "q_value",
        "p_value",
    ]
    with open(kma_result_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            record = {}
            for k, v in row.items():
                key = k.lower().replace("#", "")
                if key in int_fields:
                    try:
                        record[key] = int(v.strip())
                    except ValueError as e:
                        record[key] = None
                elif key in float_fields:
                    try:
                        record[key] = float(v.strip())
                    except ValueError as e:
                        record[key] = None
                else:
                    record[key] = v.strip()
            locus_id, allele_id, allele_hash = parse_template_name(record["template"])
            record["locus_id"] = locus_id
            record["allele_id"] = allele_id
            record["allele_hash"] = allele_hash
            if locus_id not in kma_result_by_locus_id:
                kma_result_by_locus_id[locus_id] = []
            kma_result_by_locus_id[locus_id].append(record)

    for locus_id, kma_results in kma_result_by_locus_id.items():
        kma_result_by_locus_id[locus_id] = sorted(kma_results, key=lambda k: k["score"], reverse=True)

    return kma_result_by_locus_id


BLAST_OUTFMT_FIELDS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "qlen", "slen", "nident", "bitscore"]
BLAST_OUTFMT = "6 " + " ".join(BLAST_OUTFMT_FIELDS)


def parse_blast_result(blast_result_file, min_overlap_fraction=0.5):
    """
    Parse a blastn tabular output file (outfmt BLAST_OUTFMT) into a dict of
    lists of dicts, grouped by locus_id - the same shape parse_kma_result
    produces, so downstream allele-calling logic is aligner-agnostic.

    Unlike kma's competitive read mapping, blastn reports every hit above
    its significance threshold, including several similar catalogued
    alleles all matching the very same genomic region. Hits against the
    same locus whose query ranges substantially overlap are collapsed to
    the single best-scoring one before returning, so a locus with several
    similar catalogued alleles isn't mistaken for a multi-copy locus - only
    hits at genuinely distinct query positions (or different contigs) are
    kept as separate candidates.

    depth is always None (there is no read-depth concept for an assembly).

    :param blast_result_file: Path to a blastn tabular (outfmt 6) output file
    :type blast_result_file: str
    :param min_overlap_fraction: Query ranges overlapping by at least this
        fraction of the shorter hit's length are treated as the same
        genomic copy
    :type min_overlap_fraction: float
    :return: The blast results, indexed by locus_id
    :rtype: dict[str, list[dict]]
    """
    raw_hits_by_locus_id = {}
    with open(blast_result_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t', fieldnames=BLAST_OUTFMT_FIELDS)
        for row in reader:
            locus_id, allele_id, allele_hash = parse_template_name(row["sseqid"])
            slen = int(row["slen"])
            record = {
                "template": row["sseqid"],
                "locus_id": locus_id,
                "allele_id": allele_id,
                "allele_hash": allele_hash,
                "score": float(row["bitscore"]),
                "template_length": slen,
                "template_identity": round(int(row["nident"]) / slen * 100, 2),
                "template_coverage": round((abs(int(row["send"]) - int(row["sstart"])) + 1) / slen * 100, 2),
                "query_identity": float(row["pident"]),
                "depth": None,
                "query_seqid": row["qseqid"],
                "query_start": int(row["qstart"]),
                "query_end": int(row["qend"]),
            }
            raw_hits_by_locus_id.setdefault(locus_id, []).append(record)

    return {
        locus_id: _dedupe_overlapping_hits(hits, min_overlap_fraction)
        for locus_id, hits in raw_hits_by_locus_id.items()
    }


def _dedupe_overlapping_hits(hits, min_overlap_fraction):
    kept = []
    for hit in sorted(hits, key=lambda h: h["score"], reverse=True):
        if not any(_same_genomic_region(hit, other, min_overlap_fraction) for other in kept):
            kept.append(hit)
    return kept


def _same_genomic_region(hit_a, hit_b, min_overlap_fraction):
    if hit_a["query_seqid"] != hit_b["query_seqid"]:
        return False
    overlap = min(hit_a["query_end"], hit_b["query_end"]) - max(hit_a["query_start"], hit_b["query_start"]) + 1
    if overlap <= 0:
        return False
    shorter_length = min(
        hit_a["query_end"] - hit_a["query_start"] + 1,
        hit_b["query_end"] - hit_b["query_start"] + 1,
    )
    return overlap / shorter_length >= min_overlap_fraction


def parse_kma_mapstat(kma_mapstat_file):
    """
    Parse a kma mapstat file into a dict of lists of dicts.

    :param kma_mapstat_file: The path to the kma mapstat file
    :type kma_mapstat_file: str
    :return: The kma mapstat, indexed by locus_id. Keys of kma mapstat are: locus_id, allele_id, ref_sequence, read_count, fragment_count, map_score_sum, ref_covered_positions, ref_consensus_sum, bp_total, depth_variance, nuc_high_depth_variance, depth_max, snp_sum, insert_sum, deletion_sum, read_count_aln, fragment_count_aln
    :rtype: dict[str, list[dict]]
    """
    parsed_kma_mapstat_by_locus_id = {}
    header = [
        "ref_sequence",
        "read_count",
        "fragment_count",
        "map_score_sum",
        "ref_covered_positions",
        "ref_consensus_sum",
        "bp_total",
        "depth_variance",
        "nuc_high_depth_variance",
        "depth_max",
        "snp_sum",
        "insert_sum",
        "deletion_sum",
        "read_count_aln",
        "fragment_count_aln",
    ]
    int_fields = [
        "read_count",
        "fragment_count",
        "map_score_sum",
        "ref_covered_positions",
        "ref_consensus_sum",
        "bp_total",
        "nuc_high_depth_variance",
        "depth_max",
        "snp_sum",
        "insert_sum",
        "deletion_sum",
        "read_count_aln",
        "fragment_count_aln",
    ]
    float_fields = [
        "depth_variance",
    ]
    with open(kma_mapstat_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                record = dict(zip(header, line.strip().split("\t")))
                locus_id, allele_id, allele_hash = parse_template_name(record['ref_sequence'])
                record['locus_id'] = locus_id
                record['allele_id'] = allele_id
                record['allele_hash'] = allele_hash
                for field in int_fields:
                    try:
                        record[field] = int(record[field])
                    except ValueError as e:
                        record[field] = None
                for field in float_fields:
                    try:
                        record[field] = float(record[field])
                    except ValueError as e:
                        record[field] = None
                if locus_id not in parsed_kma_mapstat_by_locus_id:
                    parsed_kma_mapstat_by_locus_id[locus_id] = []
                parsed_kma_mapstat_by_locus_id[locus_id].append(record)

    for locus_id, kma_mapstat in parsed_kma_mapstat_by_locus_id.items():
        parsed_kma_mapstat_by_locus_id[locus_id] = sorted(kma_mapstat, key=lambda k: k["map_score_sum"], reverse=True)

    return parsed_kma_mapstat_by_locus_id


def parse_kma_aln(kma_aln_file):
    """
    Parse a kma aln file into a dict of lists of dicts.
    """
    alignments_by_template_id = {}
    with open(kma_aln_file, 'r') as f:
        alignment = {}
        template_id = None
        template_seq = ""
        query_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                if template_id:
                    alignment[template_id] = {
                        'template': template_seq,
                        'query': query_seq,
                    }
                    template_seq = ""
                    query_seq = ""
                template_id = line.split(" ")[1]
            else:
                if line.startswith('template'):
                    template_seq_line = line.split(":")[1].strip()
                    template_seq += template_seq_line
                elif line.startswith('query'):
                    query_seq_line = line.split(":")[1].strip()
                    query_seq += query_seq_line

        alignment[template_id] = {
            'template': template_seq,
            'query': query_seq,
        }

    return alignments_by_template_id


def parse_allele_calls(allele_calls_path):
    """
    Parse an allele_calls.csv file back into records. Numeric fields are left
    as None for "no_hit" loci, where core-typer writes them out blank.
    """
    allele_calls = []
    int_fields = [
        'score',
        'template_length',
    ]
    float_fields = [
        'query_identity',
        'template_identity',
        'template_coverage',
        'depth',
    ]

    with open(allele_calls_path, 'r') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            for field in float_fields:
                if row[field] == '':
                    row[field] = None
                    continue
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    logger.error(f"Error parsing value: {row[field]} as float.")
                    exit(-1)
            for field in int_fields:
                if row[field] == '':
                    row[field] = None
                    continue
                try:
                    row[field] = int(row[field])
                except ValueError as e:
                    logger.error(f"Error parsing value: {row[field]} as int.")
                    exit(-1)

            allele_calls.append(row)

    return allele_calls
                

def parse_locus_names(locus_names_path):
    """
    Parse the kma index .names file to get the names of all loci
    in the order that they appear in the file.
    """
    locus_names = []
    with open(locus_names_path, 'r') as f:
        for line in f:
            locus_name, _allele_id, _allele_hash = parse_template_name(line.strip())

            if locus_name not in locus_names:
                locus_names.append(locus_name)

    return locus_names


def parse_locus_names_from_fasta(scheme_fasta_path):
    """
    Get the names of all loci in a scheme, in the order they appear, by
    reading the combined scheme FASTA's headers directly. Used for blastn
    (assembly) runs: unlike kma_index, makeblastdb doesn't produce a
    reliable, independently-readable listing of all templates (blastdbcmd
    -entry all is unreliable across blast+ versions for non-accession-style
    headers), so the original scheme FASTA is the one aligner-agnostic
    source for this.

    :param scheme_fasta_path: Path to the combined scheme FASTA (the same file passed to makeblastdb -in)
    :type scheme_fasta_path: str
    :return: Locus names, in the order they first appear
    :rtype: list[str]
    """
    locus_names = []
    with open(scheme_fasta_path, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                continue
            locus_name, _allele_id, _allele_hash = parse_template_name(line[1:].strip())

            if locus_name not in locus_names:
                locus_names.append(locus_name)

    return locus_names


def parse_kma_consensus_fasta(kma_fsa_file):
    """
    Parse a kma "-ef" consensus fasta file (kma-out.fsa) into a dict of
    template name -> consensus sequence. This is the sample's own observed
    sequence for the best-hit template at each locus, used to extract
    candidate novel alleles.

    :param kma_fsa_file: The path to the kma consensus fasta file
    :type kma_fsa_file: str
    :return: Consensus sequence, indexed by template name (e.g. "locusA_1")
    :rtype: dict[str, str]
    """
    consensus_by_template = {}
    template = None
    seq_lines = []
    with open(kma_fsa_file, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if template is not None:
                    consensus_by_template[template] = "".join(seq_lines)
                template = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if template is not None:
            consensus_by_template[template] = "".join(seq_lines)

    return consensus_by_template
            
