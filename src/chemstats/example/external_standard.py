#import chemstats
from main import ExternalStandardAnalysis


#Load 
standard_concentration = [0, 
     5, 
     10, 
     15, 
     20, 
     30, 
     35, 
     50]

signal = [
    60,        # baseline
    125000,    # ~5% con rumore
    248000,    # ~10
    372500,    # ~15
    495000,    # ~20
    742000,    # ~30
    865000,    # ~35
    1245000    # ~50
]

#
params = {
    "standards_conc": standard_concentration,
    "standards_value": signal,
    "y0": 6963644,
    "k": 3,
    "DF": 1,
    "conf": 0.95,
    "report": True,
    "report_name": "test.txt",
}

ExternalStandardAnalysis(params).run()
