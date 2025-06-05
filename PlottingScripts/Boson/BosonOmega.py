from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


mpl.rcParams.update({
    "text.usetex":       True,       # use LaTeX to render all text
    "font.family":       "serif",    # match your document font (e.g. Computer Modern)
    "font.serif":        [],         # empty means use default LaTeX serif
    "font.size":         13,         # base font size in points
    "axes.labelsize":    13,         # axis labels
    "axes.titlesize":    13,         # axis titles
    "xtick.labelsize":   13,          # tick labels
    "ytick.labelsize":   13,
    "legend.fontsize":   10,
    "legend.frameon":    True,
    "figure.dpi":        300,        # high-quality export
})

save_dir = Path('Results/Figures/Bosons/Program3')
save_dir.mkdir(parents=True, exist_ok=True)

data_dir = Path("Results/Tables/Bosons/Program7")
fname    = data_dir / "3D_Energy_omega1.000000.cvv"

df = pd.read_csv(
    fname,
    sep=r"\s+",
    comment='E',
    names=["Energy", "Omega"],  # force these column names
    header=None,                # since we're explicitly naming columns
    skip_blank_lines=True
)

df = df.astype({"Energy": float, "Omega": float})

omegas = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
analyt = 4.5 * omegas


c = sns.color_palette("colorblind")
plt.style.use("ggplot")



fig, ax = plt.subplots(figsize=(6.5, 4.8))
ax.plot(df['Omega'], df['Energy'], color=c[0], marker='o',
    label=f"Metropolis, variational energy $E_T$, 3D, $N=3$")
ax.plot(df['Omega'], analyt, 'D', mfc='white', mec=c[1], ms=6,
    label=f"analytical GS energy $E_0$, 3D, $N=3$")

ax.set_xlabel(r"Oscillator frequency $\omega$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "E_omega.pdf")



data_dir = Path("Results/Tables/Bosons/Program7")
fname    = data_dir / "3D_Energy_omega_MH1.000000.cvv"

df = pd.read_csv(
    fname,
    sep=r"\s+",
    comment='E',
    names=["Energy", "Omega"],  # force these column names
    header=None,                # since we're explicitly naming columns
    skip_blank_lines=True
)

df = df.astype({"Energy": float, "Omega": float})

omegas = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
analyt = 4.5 * omegas



fig, ax = plt.subplots(figsize=(6.5, 4.8))
ax.plot(df['Omega'], df['Energy'], color=c[2], marker='o',
    label=f"MH, variational energy $E_T$, 3D, $N=3$")
ax.plot(df['Omega'], analyt, 'D', mfc='white', mec=c[3], ms=6,
    label=f"analytical GS energy $E_0$, 3D, $N=3$")

ax.set_xlabel(r"Oscillator frequency $\omega$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "E_omega_MH.pdf")