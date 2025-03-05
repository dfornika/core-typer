import csv
import logging

from . import parsers

logger = logging.getLogger(__name__)

def choose_best_allele(kma_results, min_identity=100.0, min_coverage=100.0):
    """
    Given a list of kma results for a specific locus, choose the best allele.
    Best allele is chosen based on score. If multiple alleles have the same
    score, the first allele is chosen. If no alleles meet the minimum identity
    or coverage, the allele id is set to "-".
    
    :param kma_results: The kma results for a specific locus
    :type kma_results: list[dict]
    :param min_identity: The minimum identity required to call an allele
    :type min_identity: float
    :param min_coverage: The minimum coverage required to call an allele
    :type min_coverage: float
    :return: The best allele
    :rtype: dict
    """
    best_allele = None
    for kma_result in kma_results:
        if kma_result["template_identity"] < min_identity:
            kma_result["allele_id"] = "-"
        if kma_result["template_coverage"] < min_coverage:
            kma_result["allele_id"] = "-"
        if best_allele is None:
            best_allele = kma_result
        else:
            if kma_result["score"] > best_allele["score"]:
                best_allele = kma_result

    return best_allele


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
        "allele_id",
        "score",
        "template_length",
        "percent_identity",
        "percent_coverage",
        "depth",
    ]
    with open(allele_calls_file, 'w') as f:
        writer = csv.DictWriter(f, delimiter=',', fieldnames=output_fieldnames)
        writer.writeheader()
        for allele_call in allele_calls:
            output_record = {}
            for fieldname in output_fieldnames:
                if fieldname in allele_call:
                    output_record[fieldname] = allele_call[fieldname]
                elif fieldname == "percent_identity":
                    output_record[fieldname] = allele_call["template_identity"]
                elif fieldname == "percent_coverage":
                    output_record[fieldname] = allele_call["template_coverage"]
            writer.writerow(output_record)


def write_allele_profile(allele_calls, scheme_path, allele_profile_path):
    """
    """

    names_file = f"{scheme_path}.name"
    locus_ids = parsers.parse_locus_names(names_file)
    allele_calls_by_locus_id = {}
    for allele_call in allele_calls:
        locus_id = allele_call['locus_id']
        if locus_id not in allele_calls_by_locus_id:
            allele_calls_by_locus_id[locus_id] = allele_call
        else:
            raise ValueError(f"Duplicate allele call found for locus {locus_id}")

    with open(allele_profile_path, 'w') as f:
        f.write(','.join(locus_ids) + '\n')
        for locus_id in locus_ids:
            if locus_id in allele_calls_by_locus_id:
                allele_id = allele_calls_by_locus_id[locus_id]['allele_id']
            else:
                allele_id = '-'
            f.write(allele_id)
            if locus_id != locus_ids[-1]:
                f.write(',')
        f.write('\n')
