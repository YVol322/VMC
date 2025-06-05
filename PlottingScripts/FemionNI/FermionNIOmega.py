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

save_dir = Path('Results/Figures/Fermions/NonInteracting/Program3')
save_dir.mkdir(parents=True, exist_ok=True)

data_dir = Path("Results/Tables/Fermions/NonInteracting/Program7")
fname1    = data_dir / "N=2.cvv"
fname2    = data_dir / "N=6.cvv"
fname3    = data_dir / "N=12.cvv"

df2 = pd.read_csv(
    fname1,
    sep=r"\s+",
    comment='E',
    names=["Energy", "Omega"],  # force these column names
    header=None,                # since we're explicitly naming columns
    skip_blank_lines=True
)

df2 = df2.astype({"Energy": float, "Omega": float})

df6 = pd.read_csv(
    fname2,
    sep=r"\s+",
    comment='E',
    names=["Energy", "Omega"],  # force these column names
    header=None,                # since we're explicitly naming columns
    skip_blank_lines=True
)

df6 = df6.astype({"Energy": float, "Omega": float})

df12 = pd.read_csv(
    fname3,
    sep=r"\s+",
    comment='E',
    names=["Energy", "Omega"],  # force these column names
    header=None,                # since we're explicitly naming columns
    skip_blank_lines=True
)

df12 = df12.astype({"Energy": float, "Omega": float})

omegas = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
analyt2 = 2 * omegas
analyt6 = 10 * omegas
analyt12 = 28 * omegas


c = sns.color_palette("colorblind")
plt.style.use("ggplot")



fig, ax = plt.subplots(figsize=(6.5, 4.8))
ax.plot(df2['Omega'], df2['Energy'], color=c[0], marker='o',
    label=r"Metropolis, variational energy $E_T$, $N=2$")
ax.plot(df2['Omega'], analyt2, 'D', mfc='white', mec=c[1], ms=6,
    label=r"Analytycal GS energy $E_0=2$, $N=2$")
ax.plot(df6['Omega'], df6['Energy'], color=c[2], marker='o',
    label=r"Metropolis, variational energy $E_T$, $N=6$")
ax.plot(df6['Omega'], analyt6, 'D', mfc='white', mec=c[3], ms=6,
    label=r"Analytycal GS energy $E_0=10$, $N=6$")
ax.plot(df12['Omega'], df12['Energy'], color=c[4], marker='o',
    label=r"Metropolis, variational energy $E_T$, $N=12$")
ax.plot(df12['Omega'], analyt12, 'D', mfc='white', mec=c[7], ms=6,
    label=r"Analytycal GS energy $E_0=28$, $N=12$")

ax.set_xlabel(r"Oscillator frequency $\omega$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 10)
#plt.show()
plt.savefig(save_dir / "E_omega.pdf")





data_dir = Path("Results/Tables/Fermions/NonInteracting/Program7")
fname1    = data_dir / "MH_N=2.cvv"
fname2    = data_dir / "MH_N=6.cvv"
fname3    = data_dir / "MH_N=12.cvv"

df2 = pd.read_csv(
    fname1,
    sep=r"\s+",
    comment='E',
    names=["Energy", "Omega"],  # force these column names
    header=None,                # since we're explicitly naming columns
    skip_blank_lines=True
)

df2 = df2.astype({"Energy": float, "Omega": float})

df6 = pd.read_csv(
    fname2,
    sep=r"\s+",
    comment='E',
    names=["Energy", "Omega"],  # force these column names
    header=None,                # since we're explicitly naming columns
    skip_blank_lines=True
)

df6 = df6.astype({"Energy": float, "Omega": float})

df12 = pd.read_csv(
    fname3,
    sep=r"\s+",
    comment='E',
    names=["Energy", "Omega"],  # force these column names
    header=None,                # since we're explicitly naming columns
    skip_blank_lines=True
)

df12 = df12.astype({"Energy": float, "Omega": float})

omegas = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
analyt2 = 2 * omegas
analyt6 = 10 * omegas
analyt12 = 28 * omegas


c = sns.color_palette("colorblind")
plt.style.use("ggplot")



fig, ax = plt.subplots(figsize=(6.5, 4.8))
ax.plot(df2['Omega'], df2['Energy'], color=c[7], marker='o',
    label=r"MH, variational energy $E_T$, $N=2$")
ax.plot(df2['Omega'], analyt2, 'D', mfc='white', mec=c[6], ms=6,
    label=r"Analytycal GS energy $E_0=2$, $N=2$")
ax.plot(df6['Omega'], df6['Energy'], color=c[5], marker='o',
    label=r"MH, variational energy $E_T$, $N=6$")
ax.plot(df6['Omega'], analyt6, 'D', mfc='white', mec=c[4], ms=6,
    label=r"Analytycal GS energy $E_0=10$, $N=6$")
ax.plot(df12['Omega'], df12['Energy'], color=c[3], marker='o',
    label=r"MH, variational energy $E_T$, $N=12$")
ax.plot(df12['Omega'], analyt12, 'D', mfc='white', mec=c[2], ms=6,
    label=r"Analytycal GS energy $E_0=28$, $N=12$")

ax.set_xlabel(r"Oscillator frequency $\omega$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 10)
#plt.show()
plt.savefig(save_dir / "E_omega_MH.pdf")