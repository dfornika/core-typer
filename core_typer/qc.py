import csv
import statistics

def calculate_qc_stats(allele_calls):
    """
    """
    qc_stats = {}
    depths = []
    num_loci = 0
    num_called_alleles = 0
    for allele_call in allele_calls:
        num_loci += 1
        if allele_call['allele_id'] != '-':
            num_called_alleles += 1
        depths.append(allele_call['depth'])

    mean_depth = round(sum(depths) / len(depths), 3)
    stdev_depth = round(statistics.stdev(depths), 3)
    percent_called = round(num_called_alleles / num_loci * 100, 3)

    qc_stats['mean_depth'] = mean_depth
    qc_stats['stdev_depth'] = stdev_depth
    qc_stats['percent_called'] = percent_called

    return qc_stats

    
    
def write_qc_stats(qc_stats, qc_stats_file):
    """
    Write QC stats to a file.
    """
    qc_fieldnames = [
        'mean_depth',
        'stdev_depth',
        'percent_called',
    ]
    with open(qc_stats_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=qc_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
        writer.writeheader()
        writer.writerow(qc_stats)
