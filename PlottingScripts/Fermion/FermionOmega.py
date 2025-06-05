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

c = sns.color_palette("colorblind")
plt.style.use("ggplot")

save_dir = Path('Results/Figures/Fermions/Interacting/Program5')
save_dir.mkdir(parents=True, exist_ok=True)

data_dir = Path("Results/Tables/Fermions/Interacting/Program7")
fname1    = data_dir / "OMP_omega=0.500000_MH_PJ_eta=0.1_h=1.cvv"

df2 = pd.read_csv(
    fname1,
    sep=r"\s+",
    header=0,
    names=["Iter", "Energy", "Grad", "Beta"],
)

df2 = df2.astype({"Iter": int, "Energy": float, "Grad": float, "Beta": float})



benchmark = 1.65975

fig, ax = plt.subplots(figsize=(6.5, 4.6))
ax.plot(df2['Iter'], df2['Energy'], color=c[0], marker='o',
    label=r"Variational energy $E_T$, $\omega = 0.5$")
ax.axhline(
    y=benchmark,
    color="k",                # black line
    linestyle="--",           # dashed
    linewidth=1.5,
    label=r"DMC benchmark $E_0 = 1.65975$, $\omega = 0.5$",
)


ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
plt.legend(fontsize = 12)
plt.savefig(save_dir / "E2_omega05.pdf")


data_dir2 = Path("Results/Tables/Fermions/Interacting/Program8")
fname2    = data_dir2 / "OMP_omega=0.280000_MH_PJ_eta=0.1_h=0.1.cvv"

df6 = pd.read_csv(
    fname2,
    sep=r"\s+",
    header=0,
    names=["Iter", "Energy", "Grad"],
)

df6 = df6.astype({"Iter": int, "Energy": float, "Grad": float})


fname3    = data_dir2 / "OMP_omega=0.500000_MH_PJ_eta=0.1_h=0.1.cvv"

df66 = pd.read_csv(
    fname3,
    sep=r"\s+",
    header=0,
    names=["Iter", "Energy", "Grad"],
)

df66 = df66.astype({"Iter": int, "Energy": float, "Grad": float})



benchmark6 = 7.6001
benchmark66 = 11.7888

fig, ax = plt.subplots(figsize=(6.5, 4.6))
ax.plot(df6['Iter'], df6['Energy'], color=c[0], marker='o',
    label=r"Variational energy $E_T$, $\omega = 0.28$")
ax.axhline(
    y=benchmark6,
    color="k",                # black line
    linestyle="--",           # dashed
    linewidth=1.5,
    label=r"DMC benchmark $E_0 = 7.6001$, $\omega = 0.28$ ",
)


ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
plt.legend(fontsize = 12)
plt.savefig(save_dir / "E6_omega028.pdf")



fig, ax = plt.subplots(figsize=(7, 4.6))
ax.plot(df66['Iter'], df66['Energy'], color=c[1], marker='o',
    label=r"Variational energy $E_T$, $\omega = 0.5$")
ax.axhline(
    y=benchmark66,
    color="k",                # black line
    linestyle="--",           # dashed
    linewidth=1.5,
    label=r"DMC benchmark $E_0 = 11.7888$, $\omega = 0.5$ ",
)


ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
plt.legend(fontsize = 12)
plt.savefig(save_dir / "E6_omega05.pdf")
#plt.show()





data_dir3 = Path("Results/Tables/Fermions/Interacting/Program9")
fname4    = data_dir3 / "OMP_omega=0.280000_MH_PJ_eta=0.01_h=0.1.cvv"

df12 = pd.read_csv(
    fname4,
    sep=r"\s+",
    header=0,
    names=["Iter", "Energy", "Grad", "Beta"],
)

df12 = df12.astype({"Iter": int, "Energy": float, "Grad": float, "Beta": float})

benchmark12 = 25.6356

fig, ax = plt.subplots(figsize=(6.5, 4.6))
ax.plot(df12['Iter'], df12['Energy'], color=c[0], marker='o',
    label=r"Variational energy $E_T$, $\omega = 0.28$")
ax.axhline(
    y=benchmark12,
    color="k",                # black line
    linestyle="--",           # dashed
    linewidth=1.5,
    label=r"DMC benchmark $E_0 = 25.6356$, $\omega = 0.28$ ",
)


ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
plt.legend(fontsize = 12)
plt.savefig(save_dir / "E12_omega028.pdf")
#plt.show()



fname5    = data_dir3 / "OMP_omega=0.500000_MH_PJ_eta=0.01_h=0.1.cvv"

df122 = pd.read_csv(
    fname5,
    sep=r"\s+",
    header=0,
    names=["Iter", "Energy", "Grad", "Beta"],
)

df122 = df122.astype({"Iter": int, "Energy": float, "Grad": float, "Beta": float})


benchmark12 = 39.159

fig, ax = plt.subplots(figsize=(6.5, 4.6))
ax.plot(df122['Iter'], df122['Energy'], color=c[1], marker='o',
    label=r"Variational energy $E_T$, $\omega = 0.5$")
ax.axhline(
    y=benchmark12,
    color="k",                # black line
    linestyle="--",           # dashed
    linewidth=1.5,
    label=r"DMC benchmark $E_0 = 39.159$, $\omega = 0.5$ ",
)


ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
plt.legend(fontsize = 12)
plt.savefig(save_dir / "E12_omega05.pdf")
#plt.show()