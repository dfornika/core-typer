import csv

def parse_kma_result(kma_result_file):
    """
    Parse a kma result file into a dict of lists of dicts.

    :param kma_result_file: The path to the kma result file
    :type kma_result_file: str
    :return: The kma results, indexed by locus_id. Keys of kma results are: locus_id, allele_id, score, template_length, template_identity, template_coverage, query_identity, query_coverage, depth, q_value, p_value
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
                locus_id = record["template"].split("_")[0]
                allele_id = record["template"].split("_")[1]
                record["locus_id"] = locus_id
                record["allele_id"] = allele_id
                if locus_id not in kma_result_by_locus_id:
                    kma_result_by_locus_id[locus_id] = []
                kma_result_by_locus_id[locus_id].append(record)

    for locus_id, kma_results in kma_result_by_locus_id.items():
        kma_result_by_locus_id[locus_id] = sorted(kma_results, key=lambda k: k["score"], reverse=True)

    return kma_result_by_locus_id


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
                locus_id = record['ref_sequence'].split("_")[0]
                allele_id = record['ref_sequence'].split("_")[1]
                record['locus_id'] = locus_id
                record['allele_id'] = allele_id
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
        for line in f:
            if line.startswith("#"):
                template_id = line.strip().split(" ")[1]
            else:
                if line.startswith('template'):
                    sequence = line.strip().split(" ")[1]
