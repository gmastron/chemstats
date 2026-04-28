import matplotlib.pyplot as plt
import numpy as np
from .ExternalStandard import ExternalStandard


class ExternalStandardAnalysis:

    def __init__(self, params):
        self.model = ExternalStandard(params)

    def run(self):
        self.model.fit()
        self.model.calc_syx()
        self.model.plot_regression_line()
        self.model.interpolate()

        if self.model.conf is not None:
            self.model.confidence_evaluation()
            self.model.print_report()
            
            