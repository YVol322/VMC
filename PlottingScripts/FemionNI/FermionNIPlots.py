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


save_dir = Path('Results/Figures/Fermions/NonInteracting')
save_dir.mkdir(parents=True, exist_ok=True)

bases = {
    'Serial': Path('Results/Tables/Fermions/NonInteracting/Program1'),
    'MPI':    Path('Results/Tables/Fermions/NonInteracting/Program2'),
    'OMP':    Path('Results/Tables/Fermions/NonInteracting/Program3'),
}

parts = ['2', '6', '12']
cols  = ['Energy', 'Time']

data = {}

for backend, base_path in bases.items():
    # decide the filename prefix
    prefix = '' if backend == 'Serial' else f"{backend}_"

    for part in parts:
        fname = f"{prefix}N={part}.cvv"
        path  = base_path / fname

        if not path.exists():
            raise FileNotFoundError(f"Missing file: {path}")

        df = pd.read_csv(path, sep='\t')[cols]
        # store under (backend, part)
        data[(backend, part)] = df


E_analyt = [2, 10, 28]
E_Metro = np.array([data['Serial','2'].Energy[0], data['Serial','6'].Energy[0], data['Serial','12'].Energy[0]])
E_Num = np.array([data['Serial','2'].Energy[1], data['Serial','6'].Energy[1], data['Serial','12'].Energy[1]])
E_MH = np.array([data['Serial','2'].Energy[2], data['Serial','6'].Energy[2], data['Serial','12'].Energy[2]])
E_NumMH = np.array([data['Serial','2'].Energy[3], data['Serial','6'].Energy[3], data['Serial','12'].Energy[3]])


TMetro2 = data['Serial','2'].Time[0]
TMetroNum2 = data['Serial','2'].Time[1]
TMH2 = data['Serial','2'].Time[2]
TMHNum2 = data['Serial','2'].Time[3]

TMetro6 = data['Serial','6'].Time[0]
TMetroNum6 = data['Serial','6'].Time[1]
TMH6 = data['Serial','6'].Time[2]
TMHNum6 = data['Serial','6'].Time[3]

TMetro12 = data['Serial','12'].Time[0]
TMetroNum12 = data['Serial','12'].Time[1]
TMH12 = data['Serial','12'].Time[2]
TMHNum12 = data['Serial','12'].Time[3]

TOMPMetro2 = data['OMP','2'].Time[0]
TOMPMetroNum2 = data['OMP','2'].Time[1]
TOMPMH2 = data['OMP','2'].Time[2]
TOMPMHNum2 = data['OMP','2'].Time[3]

TOMPMetro6 = data['OMP','6'].Time[0]
TOMPMetroNum6 = data['OMP','6'].Time[1]
TOMPMH6 = data['OMP','6'].Time[2]
TOMPMHNum6 = data['OMP','6'].Time[3]

TOMPMetro12 = data['OMP','12'].Time[0]
TOMPMetroNum12 = data['OMP','12'].Time[1]
TOMPMH12 = data['OMP','12'].Time[2]
TOMPMHNum12 = data['OMP','12'].Time[3]


TMPIMetro2 = data['MPI','2'].Time[0]
TMPIMetroNum2 = data['MPI','2'].Time[1]
TMPIMH2 = data['MPI','2'].Time[2]
TMPIMHNum2 = data['MPI','2'].Time[3]

TMPIMetro6 = data['MPI','6'].Time[0]
TMPIMetroNum6 = data['MPI','6'].Time[1]
TMPIMH6 = data['MPI','6'].Time[2]
TMPIMHNum6 = data['MPI','6'].Time[3]

TMPIMetro12 = data['MPI','12'].Time[0]
TMPIMetroNum12 = data['MPI','12'].Time[1]
TMPIMH12 = data['MPI','12'].Time[2]
TMPIMHNum12 = data['MPI','12'].Time[3]


TMetro = np.array([TMetro2, TMetro6, TMetro12])
TMetroNum = np.array([TMetroNum2, TMetroNum6, TMetroNum12])
TMH = np.array([TMH2, TMH6, TMH12])
TMHNum = np.array([TMHNum2, TMHNum6, TMHNum12])

TMetroOMP = [TOMPMetro2, TOMPMetro6, TOMPMetro12]
TMetroNumOMP = [TOMPMetroNum2, TOMPMetroNum6, TOMPMetroNum12]
TMHOMP = [TOMPMH2, TOMPMH6, TOMPMH12]
TMHNumOMP = [TOMPMHNum2, TOMPMHNum6, TOMPMHNum12]

TMetroMPI = [TMPIMetro2, TMPIMetro6, TMPIMetro12]
TMetroNumMPI = [TMPIMetroNum2, TMPIMetroNum6, TMPIMetroNum12]
TMHMPI = [TMPIMH2, TMPIMH6, TMPIMH12]
TMHNumMPI = [TMPIMHNum2, TMPIMHNum6, TMPIMHNum12]


nparts = np.array([2, 6, 12])


colors = sns.color_palette("colorblind")
plt.style.use("ggplot")

c1 = colors[0]
c2 = colors[1]
c3 = colors[2]
c4 = colors[3]
c5 = colors[4]
c6 = colors[5]
c7 = colors[6]
c8 = colors[7]

fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, E_Metro - E_analyt, color=c1, marker='o', label=r"Metropolis, $c=0$")
ax.plot(nparts, E_Num - E_analyt + 1, color=c2, marker='o', label=r"Metropolis numerical, $c=1$")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Energy difference $\Delta E + c$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "E_Metro.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, E_MH - E_analyt + 2, color=c3, marker='o', label=r"MH, $c=2$")
ax.plot(nparts, E_NumMH - E_analyt + 3, color=c4, marker='o', label=r"MH numerical, $c=3$")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Energy difference $\Delta E + c$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "E_MH.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, TMetro, color=c1, marker='o', label=r"Metropolis")
ax.plot(nparts, TMetroNum, color=c2, marker='o', label=r"Metropolis numerical")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "T_Metro.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, TMH, color=c3, marker='o', label=r"MH")
ax.plot(nparts, TMHNum, color=c4, marker='o', label=r"MH numerical")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "T_MH.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, TMetro, color=c1, marker='o', label=r"Serial Metropolis")
ax.plot(nparts, TMetroOMP, color=c2, linestyle=':', marker='o', label=r"OMP Metropolis")
ax.plot(nparts, TMetroMPI, color=c3, marker='o', label=r"MPI Metropolis")

ax.set_xlabel(r"Number of particles $N$")
ax.set_ylabel(r"Time $t$, $[s]$")
ax.legend(fontsize = 12)
plt.savefig(save_dir / "T_Metro.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, TMetroNum, color=c4, marker='o', label=r"Serial Metropolis numerical")
ax.plot(nparts, TMetroNumOMP, color=c5, linestyle=':', marker='o', label=r"OMP Metropolis numerical")
ax.plot(nparts, TMetroNumMPI, color=c6, marker='o', label=r"MPI Metropolis numerical")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "T_MetroNum.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, TMH, color=c1, marker='o', label=r"Serial MH")
ax.plot(nparts, TMHOMP, color=c2, linestyle=':', marker='o', label=r"OMP MH")
ax.plot(nparts, TMHMPI, color=c3, marker='o', label=r"MPI, MH")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "T_MH.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, TMHNum, color=c4, marker='o', label=r"Serial MH numerical")
ax.plot(nparts, TMHNumOMP, color=c5, linestyle=':', marker='o', label=r"OMP MH numerical")
ax.plot(nparts, TMHNumMPI, color=c6, marker='o', label=r"MPI MH numerical")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "T_MHNum.pdf")