"""Current native-scan claims from the complete 347-cache CapScan comparison."""
from pathlib import Path

import pandas as pd

from text_statistics.supplementary_methods_6_7 import stats
from tools.check_updated_analyses import validate_roster

REPO = Path(__file__).resolve().parents[1]


def numbers():
    acquisitions = pd.read_csv(REPO / 'source_data/v42_updates/crossfit_per_acquisition.tsv',
                               sep='\t', float_precision='round_trip')
    validate_roster(acquisitions)
    values = stats.stat_crossfit(REPO / 'source_data/supplementary_methods_6_7',
                                acquisitions=acquisitions)
    return {key: stats.recomputed(values, key, printed)
            for line, quote, printed, key, function, source, status in stats.NUMBER_SPEC
            if 'crossfit_per_acquisition.tsv' in source}
