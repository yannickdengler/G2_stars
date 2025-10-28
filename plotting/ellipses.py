import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

class measurement:
    def __init__(self, name, R, M, R_err, M_err, text_R, text_M):
        # these are instance attributes
        self.name = name
        self.R = R
        self.M = M
        self.R_err = R_err
        self.M_err = M_err
        self.text_R = text_R
        self.text_M = text_M
    def draw_ellips(self, ax, color = "red"):
        ellipse = Ellipse((self.R, self.M), width=2*self.R_err, height=2*self.M_err, angle=0,
                        edgecolor=None, facecolor=color, alpha=0.5,zorder = 11)
        ax.add_patch(ellipse)
    def write_text(self, ax):
        props = dict(facecolor='white', alpha=0.5, edgecolor="None")#, boxstyle="round")
        ax.text(x = self.text_R, y = self.text_M, s = self.name, bbox=props, zorder = 10)

# HESS = measurement("HESS J1731-347"     , 10.40, 0.770, 0.82, 0.188, 10.8, 0.660)
# PSR_light = measurement("PSR J0030+0451", 13.11, 1.370, 1.30, 0.170, 13.11, 1.370)
# PSR_mid = measurement("PSR J0437-4715"  , 11.36, 1.418, 0.79, 0.037, 11.36, 1.45)
# PSR_heavy = measurement("PSR J0740+6620", 12.49, 2.073, 1.08, 0.069, 11.8, 2.14)
# GW170817 = measurement("GW170817"       , 10.75, 1.360, 1.11, 0.215, 8.4, 1.245)
HESS = measurement(r"$\text{HESS J1731-347}$"     , 10.40, 0.770, 0.82, 0.185, 8.3, 0.97)
# PSR_light = measurement("PSR J0030+0451", 13.11, 1.370, 1.30, 0.170, 13.11, 1.370)
PSR_mid = measurement(r"$\text{PSR J0437-4715}$"  , 11.36, 1.418, 0.79, 0.037, 8.5, 1.61)
PSR_heavy = measurement(r"$\text{PSR J0740+6620}$", 12.49, 2.073, 1.08, 0.069, 10.1, 2.17)
# GW170817 = measurement(r"$\text{GW170817}$"       , 11.90, 1.385, 1.40, 0.215, 8.95, 1.2)
# GW170817 = measurement(r"$\text{GW170817}$"       , 11.90, 1.385, 1.40, 0.215, 12.7, 1.12)
GW170817 = measurement(r"$\text{GW170817}$"       , 11.90, 1.385, 1.40, 0.215, 9, 1.2)

c_14 = "olivedrab"
c_10_res = "royalblue"
c_10_non_res = "orangered"
def draw_errors(ax):
    HESS.draw_ellips(ax, "royalblue")
    PSR_heavy.draw_ellips(ax, "mediumorchid")
    # GW170817.draw_ellips(ax, "orangered")
    GW170817.draw_ellips(ax, "goldenrod")
    PSR_mid.draw_ellips(ax, "olivedrab")

def draw_error_text(ax):
    ax.annotate("", [11.2,1.4], [10,1.6],
            arrowprops=dict(arrowstyle="->"),zorder=12)
    HESS.write_text(ax)
    PSR_heavy.write_text(ax)
    GW170817.write_text(ax)
    PSR_mid.write_text(ax)

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize = (8,6))
    draw_errors(ax)

    ax.set_xlim(0, 15)
    ax.set_ylim(0, 2.5)
    plt.show()
