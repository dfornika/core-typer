import csv
import hashlib
import logging

logger = logging.getLogger(__name__)

CLEAN_BASES = {"A", "C", "G", "T"}


def choose_best_allele(kma_results, min_identity=100.0, min_coverage=100.0):
    """
    Given a list of kma results for a specific locus, choose the best allele.
    Best allele is chosen based on score. If multiple alleles have the same
    score, the first allele is chosen.

    The identity check uses query_identity (mismatches within the aligned
    region), not template_identity (which conflates mismatches with missing
    coverage into one number). This matches CGE cgMLSTFinder's approach and
    keeps "diverged from the closest allele" and "gene only partially
    assembled" distinguishable via query_identity vs. template_coverage.

    The returned record is always annotated with:
      - closest_allele_id / closest_allele_hash: the allele id (and, for
        schemes prepared with allele hashing, md5 hash) of the best-scoring
        hit, kept even when it doesn't meet the thresholds, so the nearest
        known allele isn't lost
      - call_status: "called" if the best hit meets min_identity and
        min_coverage, otherwise "divergent"
      - allele_id: the called allele id if call_status is "called", else "-"
      - num_hits: total number of candidate templates kma reported for
        this locus. A genuine second (or third...) copy of a locus, if
        divergent enough from the first, can make kma report more than one
        hit here - evidence of a possible multi-copy locus that this
        function would otherwise silently discard by keeping only the best
        hit. See find_possible_multicopy_loci().

    :param kma_results: The kma results for a specific locus
    :type kma_results: list[dict]
    :param min_identity: The minimum identity required to call an allele
    :type min_identity: float
    :param min_coverage: The minimum coverage required to call an allele
    :type min_coverage: float
    :return: The best allele
    :rtype: dict
    """
    best_allele = kma_results[0]
    for kma_result in kma_results[1:]:
        if kma_result["score"] > best_allele["score"]:
            best_allele = kma_result

    best_allele["closest_allele_id"] = best_allele["allele_id"]
    best_allele["closest_allele_hash"] = best_allele.get("allele_hash")
    # The original, un-split KMA template name (e.g. "locusA_1" or, for
    # hash-prepared schemes, "locusA_1_<md5>") - needed as-is to look the
    # locus up in the kma "-ef" consensus (.fsa) output, which is keyed by
    # the exact template string, hash suffix included.
    best_allele["closest_allele_template"] = best_allele.get("template")
    best_allele["novel_allele_hash"] = None
    best_allele["num_hits"] = len(kma_results)
    if best_allele["query_identity"] >= min_identity and best_allele["template_coverage"] >= min_coverage:
        best_allele["call_status"] = "called"
    else:
        best_allele["call_status"] = "divergent"
        best_allele["allele_id"] = "-"

    return best_allele


def build_complete_allele_calls(locus_ids, parsed_kma_result, min_identity=100.0, min_coverage=100.0):
    """
    Build one allele call record per locus in the scheme, in scheme order.
    Loci with no KMA hits at all (not just below-threshold hits) get a
    placeholder record with call_status "no_hit", allele_id "-", and other
    fields set to None.

    :param locus_ids: All locus ids in the scheme, in scheme order
    :type locus_ids: list[str]
    :param parsed_kma_result: The kma results, indexed by locus_id
    :type parsed_kma_result: dict[str, list[dict]]
    :param min_identity: The minimum identity required to call an allele
    :type min_identity: float
    :param min_coverage: The minimum coverage required to call an allele
    :type min_coverage: float
    :return: One allele call record per locus, in scheme order
    :rtype: list[dict]
    """
    allele_calls = []
    for locus_id in locus_ids:
        if locus_id in parsed_kma_result:
            best_allele = choose_best_allele(parsed_kma_result[locus_id], min_identity=min_identity, min_coverage=min_coverage)
        else:
            best_allele = {
                "locus_id": locus_id,
                "allele_id": "-",
                "closest_allele_id": None,
                "closest_allele_hash": None,
                "closest_allele_template": None,
                "novel_allele_hash": None,
                "num_hits": 0,
                "call_status": "no_hit",
                "score": None,
                "template_length": None,
                "template_identity": None,
                "template_coverage": None,
                "query_identity": None,
                "depth": None,
            }
        allele_calls.append(best_allele)

    return allele_calls


def find_possible_multicopy_loci(parsed_kma_result):
    """
    For each locus where kma reported more than one candidate hit, return
    all of that locus's hits as possible multi-copy evidence - independent
    of whether any of them individually meet the calling thresholds. A
    genuine second (or third...) copy of a locus, if divergent enough from
    the first, can make kma report it as a separate hit rather than merging
    it into the first; build_complete_allele_calls/choose_best_allele only
    ever keep the single best-scoring hit per locus, so this is currently
    the only place that evidence survives.

    :param parsed_kma_result: The kma results, indexed by locus_id
    :type parsed_kma_result: dict[str, list[dict]]
    :return: One record per hit, for loci with more than one hit, sorted by locus_id then score descending
    :rtype: list[dict]
    """
    multicopy_records = []
    for locus_id in sorted(parsed_kma_result):
        kma_results = parsed_kma_result[locus_id]
        if len(kma_results) <= 1:
            continue
        for kma_result in sorted(kma_results, key=lambda r: r["score"], reverse=True):
            multicopy_records.append({
                "locus_id": locus_id,
                "allele_id": kma_result["allele_id"],
                "query_identity": kma_result["query_identity"],
                "template_coverage": kma_result["template_coverage"],
                "depth": kma_result["depth"],
                "score": kma_result["score"],
            })

    return multicopy_records


def write_possible_multicopy_loci(multicopy_records, multicopy_loci_path):
    """
    Write possible multi-copy locus evidence to a CSV file, one row per hit.

    :param multicopy_records: Records from find_possible_multicopy_loci
    :type multicopy_records: list[dict]
    :param multicopy_loci_path: The path to write the report to
    :type multicopy_loci_path: str
    :return: None
    """
    fieldnames = ["locus_id", "allele_id", "query_identity", "template_coverage", "depth", "score"]
    with open(multicopy_loci_path, 'w') as f:
        writer = csv.DictWriter(f, delimiter=',', fieldnames=fieldnames)
        writer.writeheader()
        for record in multicopy_records:
            writer.writerow(record)


def is_clean_sequence(sequence):
    """
    True if sequence contains only unambiguous bases (A/C/G/T). Used to
    decide whether an observed consensus sequence is trustworthy enough to
    treat as a candidate novel allele.
    """
    return set(sequence) <= CLEAN_BASES


def hash_sequence(sequence):
    """
    Content hash for an allele sequence, used both to embed a stable
    identifier in a prepared scheme (see scripts/prepare_cgmlst_scheme.py)
    and to fingerprint an observed novel-allele candidate for the same
    comparison. MD5 is used only as a content identifier here (not for any
    security purpose), matching the convention CGE cgMLSTFinder itself uses.
    """
    return hashlib.md5(sequence.encode("utf-8")).hexdigest()


def extract_novel_alleles(allele_calls, consensus_by_template, min_coverage=100.0):
    """
    For each "divergent" locus (best hit found, but below the identity/
    coverage thresholds) that nonetheless has full template coverage, pull
    its observed consensus sequence out of kma's "-ef" fsa output and hash
    it - a candidate for submission as a genuinely new allele. Requires full
    coverage (not just whatever --min-coverage the user configured for
    general calling) since a novel allele candidate needs to be complete to
    be submittable, and requires a clean (unambiguous) consensus, matching
    what CGE cgMLSTFinder does for its "hypothetical new allele" output.

    Mutates each qualifying allele_call in place, setting novel_allele_hash.

    :param allele_calls: One allele call record per locus, in scheme order
    :type allele_calls: list[dict]
    :param consensus_by_template: Consensus sequence, indexed by template name (from parsers.parse_kma_consensus_fasta)
    :type consensus_by_template: dict[str, str]
    :param min_coverage: Minimum template_coverage required to treat a divergent locus as a novel-allele candidate
    :type min_coverage: float
    :return: One record per novel-allele candidate found
    :rtype: list[dict]
    """
    novel_alleles = []
    for allele_call in allele_calls:
        if allele_call["call_status"] != "divergent":
            continue
        if allele_call["template_coverage"] is None or allele_call["template_coverage"] < min_coverage:
            continue

        consensus = consensus_by_template.get(allele_call["closest_allele_template"])
        if consensus is None or not is_clean_sequence(consensus):
            continue

        novel_allele_hash = hash_sequence(consensus)
        allele_call["novel_allele_hash"] = novel_allele_hash
        novel_alleles.append({
            "locus_id": allele_call["locus_id"],
            "hash": novel_allele_hash,
            "sequence": consensus,
        })

    return novel_alleles


def write_novel_alleles(novel_alleles, novel_alleles_path):
    """
    Write candidate novel allele sequences to a FASTA file, one entry per
    locus, headers carrying the locus id and its content hash.

    :param novel_alleles: Novel allele records from extract_novel_alleles
    :type novel_alleles: list[dict]
    :param novel_alleles_path: The path to write the novel alleles fasta to
    :type novel_alleles_path: str
    :return: None
    """
    with open(novel_alleles_path, 'w') as f:
        for novel_allele in novel_alleles:
            f.write(f">{novel_allele['locus_id']} novel_allele_hash={novel_allele['hash']}\n")
            f.write(novel_allele['sequence'] + "\n")


def write_allele_calls(allele_calls_file, allele_calls):
    """
    Write allele calls to a CSV file.

    :param allele_calls_file: The path to the allele calls file
    :type allele_calls_file: str
    :param allele_calls: The allele calls
    :type allele_calls: list[dict]
    :return: None
    """
    output_fieldnames = [
        "locus_id",
        "call_status",
        "allele_id",
        "closest_allele_id",
        "closest_allele_hash",
        "novel_allele_hash",
        "query_identity",
        "template_identity",
        "template_coverage",
        "depth",
        "score",
        "template_length",
        "num_hits",
    ]
    with open(allele_calls_file, 'w') as f:
        writer = csv.DictWriter(f, delimiter=',', fieldnames=output_fieldnames)
        writer.writeheader()
        for allele_call in allele_calls:
            output_record = {fieldname: allele_call.get(fieldname) for fieldname in output_fieldnames}
            writer.writerow(output_record)


def write_allele_profile(allele_calls, allele_profile_path):
    """
    Write a single-row allele profile (locus ids as header, allele ids as the row),
    in the same order as allele_calls.

    :param allele_calls: One allele call record per locus, in scheme order
    :type allele_calls: list[dict]
    :param allele_profile_path: The path to write the allele profile to
    :type allele_profile_path: str
    :return: None
    """
    with open(allele_profile_path, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow([allele_call['locus_id'] for allele_call in allele_calls])
        writer.writerow([allele_call['allele_id'] for allele_call in allele_calls])
