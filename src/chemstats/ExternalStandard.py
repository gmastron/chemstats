# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
External standard calibration model
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress, t


class ExternalStandard:

    def __init__(self, params):
        self.std_conc = np.array(params["standards_conc"])
        self.std_value = np.array(params["standards_value"])
        self.y0 = params["y0"]
        self.k = params["k"]
        self.DF = params["DF"]
        self.conf = params["conf"]
        self.report = params["report"]
        self.report_name = params["report_name"]
        
        # placeholders
        self.m = None
        self.q = None
        self.r2 = None
        self.syx = None
        self.n = None
        self.x0 = None
        self.sx0_value = None
     
    # -------------------------
    # STATIC METHODS (math)
    # -------------------------
    @staticmethod
    def linear_fit(x, y):
        return linregress(x, y)

    @staticmethod
    def sx0(syx, m, q, x, y0, k):
        x = np.array(x)

        n = len(x)
        x_mean = np.mean(x)
        Sxx = np.sum((x - x_mean) ** 2)

        x0 = (y0 - q) / m

        if Sxx == 0:
            raise ValueError("Sxx = 0 (all x values identical)")

        term = (1 / k) + (1 / n) + ((x0 - x_mean) ** 2) / Sxx
        sx0 = (syx / m) * np.sqrt(term)

        return sx0, x0

    @staticmethod
    def confidence_interval(x0, sx0, n, confidence):
        df = n - 2
        t_value = t.ppf((1 + confidence) / 2, df)
        margin = t_value * sx0
        return x0 - margin, x0 + margin, margin

    @staticmethod
    def line_definition(m, q, x):
        x_line = np.linspace(min(x), max(x), 100)
        y_line = m * x_line + q
        return x_line, y_line

    # -------------------------
    # MODEL METHODS
    # -------------------------
    def fit(self):
        result = self.linear_fit(self.std_conc, self.std_value)

        self.m = result.slope
        self.q = result.intercept
        self.r2 = result.rvalue ** 2

        return self.m, self.q, self.r2

    def calc_syx(self):
        y_hat = self.m * self.std_conc + self.q
        self.n = len(self.std_conc)

        self.syx = np.sqrt(
            np.sum((self.std_value - y_hat) ** 2) / (self.n - 2)
        )

        return self.syx

    def plot_regression_line(self):
        x_line, y_line = self.line_definition(
            self.m, self.q, self.std_conc
        )

        plt.figure(figsize=(7, 7))
        plt.plot(x_line, y_line, '--', label="Regression line")
        plt.plot(self.std_conc, self.std_value, 'o', label="Standards")

        plt.text(
            min(self.std_conc),
            max(self.std_value),
            f"y = {self.m:.5f}x + {self.q:.5f}\nR² = {self.r2:.5f}"
        )

        plt.xlabel("Concentration")
        plt.ylabel("Signal")
        plt.legend()
        
        if self.report:
            plot_name = self.report_name.replace(".txt", ".png")
            plt.savefig(plot_name, dpi=300, bbox_inches="tight")
        
        
        plt.show()

    def interpolate(self):
        self.sx0_value, self.x0 = self.sx0(
            self.syx,
            self.m,
            self.q,
            self.std_conc,
            self.y0,
            self.k
        )

        self.x0_real = self.x0 * self.DF
        self.sx0_real = self.sx0_value * self.DF
        
        if self.conf is None:
            print(f"x0 = {self.x0_real:.4g}")
            
        return self.x0_real, self.sx0_real

    def confidence_evaluation(self):
        if self.conf is None:
            print(f"x0 = {self.x0_real:.4g}")
            return

        low, high, margin = self.confidence_interval(
            self.x0_real,
            self.sx0_real,
            self.n,
            self.conf
        )
        
        self.low = low
        self.high = high
        self.margin = margin
        
        err = (margin / self.x0_real) * 100
        conf_pct = self.conf * 100
        self.err = err
        self.conf_pct = conf_pct
        print(f"x0 = {self.x0_real:.4g} ± {margin:.2g} (CI {conf_pct:.1f}%)")
        print(f"Interval: [{low:.4g}, {high:.4g}]")
        print(f"Error (%): {err:.2g}")
        
    def print_report(self):
        clean_name = self.report_name.replace(".txt", "")
        if not self.report:
            return

        with open(self.report_name, "w") as f:

            f.write(f"{('CALIBRATION REPORT: ' + clean_name):^50}\n")
            f.write("=" * 50 + "\n\n")
            f.write("\n")
            f.write(f"{'CALIBRATION STANDARDS':^50}\n")
            f.write("-" * 50 + "\n")
            f.write(f"{'Concentration':<20}{'Signal':>20}\n")
            f.write("-" * 50 + "\n")
            for conc, signal in zip(self.std_conc, self.std_value):
                f.write(f"{conc:<20.4f}{signal:>20.4f}\n")
                
            f.write("=" * 50 + "\n\n")
            f.write(f"{'Parameters':^50}\n")
            f.write("-" * 50 + "\n")
            f.write(f"{'Slope:':<15}{self.m:>20.6f}\n")
            f.write(f"{'Intercept:':<15}{self.q:>20.6f}\n")
            f.write(f"{'R2:':<15}{self.r2:>20.6f}\n")
            f.write(f"{'SD:':<15}{self.sx0_real:>20.6f}\n")
            f.write(f"{'x0:':<15}{self.x0_real:>20.6f}\n")
            f.write(f"{'Margin:':<15}{self.margin:>20.6f}\n")
            f.write(f"{'Confidence:':<15}{self.conf_pct:>19.1f}%\n\n")
            f.write(f"{'OUTPUT':^50}\n")
            f.write("-" * 50 + "\n")
            f.write(f"x0 = {self.x0_real:.4g} ± {self.margin:.2g}\n")
            f.write(f"Interval: [{self.low:.4g}, {self.high:.4g}]\n")
            f.write(f"Error (%): {self.err:.2g}\n")
