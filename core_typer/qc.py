import csv
import statistics


def calculate_qc_stats(allele_calls, possible_multicopy_loci=None):
    """
    Calculate QC stats from a complete list of allele calls (one record per
    locus in the scheme, including loci with no KMA hit at all).

    mean_depth/stdev_depth are calculated only over loci that received at
    least one KMA hit (called or not), since loci with no hit have no depth
    measurement. stdev_depth requires at least 2 such loci.

    num_possible_multicopy_loci counts the distinct loci reported by
    allele_calling.find_possible_multicopy_loci, so it matches
    possible_multicopy_loci.csv (including its depth-ratio filtering for read
    input). When possible_multicopy_loci is not supplied it falls back to
    counting loci whose best call carried more than one raw candidate hit
    (num_hits > 1), which is unfiltered.

    :param allele_calls: One allele call record per locus in the scheme
    :type allele_calls: list[dict]
    :param possible_multicopy_loci: Records from find_possible_multicopy_loci
        (one per hit); distinct locus_ids are counted
    :type possible_multicopy_loci: list[dict] or None
    :return: QC stats
    :rtype: dict
    """
    num_loci = len(allele_calls)
    num_called_alleles = sum(1 for allele_call in allele_calls if allele_call['call_status'] == 'called')
    num_divergent = sum(1 for allele_call in allele_calls if allele_call['call_status'] == 'divergent')
    num_no_hit = sum(1 for allele_call in allele_calls if allele_call['call_status'] == 'no_hit')
    if possible_multicopy_loci is None:
        num_possible_multicopy_loci = sum(1 for allele_call in allele_calls if allele_call.get('num_hits', 0) > 1)
    else:
        num_possible_multicopy_loci = len({record['locus_id'] for record in possible_multicopy_loci})
    depths = [allele_call['depth'] for allele_call in allele_calls if allele_call.get('depth') is not None]

    qc_stats = {
        'num_loci': num_loci,
        'num_called_alleles': num_called_alleles,
        'num_divergent': num_divergent,
        'num_no_hit': num_no_hit,
        'num_possible_multicopy_loci': num_possible_multicopy_loci,
        'percent_called': round(num_called_alleles / num_loci * 100, 3) if num_loci else 0.0,
        'mean_depth': round(sum(depths) / len(depths), 3) if depths else None,
        'stdev_depth': round(statistics.stdev(depths), 3) if len(depths) > 1 else None,
    }

    return qc_stats


def write_qc_stats(qc_stats, qc_stats_file):
    """
    Write QC stats to a file.
    """
    qc_fieldnames = [
        'num_loci',
        'num_called_alleles',
        'num_divergent',
        'num_no_hit',
        'num_possible_multicopy_loci',
        'percent_called',
        'mean_depth',
        'stdev_depth',
    ]
    with open(qc_stats_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=qc_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
        writer.writeheader()
        writer.writerow(qc_stats)
