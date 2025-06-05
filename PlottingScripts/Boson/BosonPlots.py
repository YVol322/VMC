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

save_dir = Path('Results/Figures/Bosons/Program1')
save_dir.mkdir(parents=True, exist_ok=True)


bases = {
    'Seq': Path('Results/Tables/Bosons/Program1'),
    'MPI': Path('Results/Tables/Bosons/Program2'),
    'OMP': Path('Results/Tables/Bosons/Program3'),
}

dims = ['1D', '2D', '3D']
schemes = {
    'Metro':        '',
    'NumMetro':     'Num_',
    'MetroHast':    'MH_',
    'NumMetroHast': 'Num_MH_',
}
cols = ['N','Energy','Time','AccRatio']

data = {}

for backend, base_path in bases.items():
    for kind, prefix in schemes.items():
        for dim in dims:
            fname = f"{dim}_{prefix}Energy_Time_AccRate.cvv"
            path  = base_path / fname

            # optional: check existence
            if not path.exists():
                raise FileNotFoundError(f"{backend} file not found: {path!r}")

            df = pd.read_csv(path, sep='\t')[cols]
            data[(backend, kind, dim)] = df
        

E_analyt = {
    '1D': np.arange(0.5, 5.1, 0.5),
    '2D': np.arange(1,   10.1, 1),
    '3D': np.arange(1.5, 15.1, 1.5),
}


colors = sns.color_palette("colorblind")
plt.style.use("ggplot")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
dims = ['1D', '2D', '3D']
for i, dim in enumerate(dims):
    df = data[('Seq', 'Metro', dim)]
    c = colors[i]
    ax.plot(df['N'], df['Energy'], color=c, marker='o',
            label=f"Metropolis variational energy $E_T$, {dim}")
    ax.plot(df['N'], E_analyt[dim], 'D', mfc='white', mec=c, ms=6,
            label=f"Metropolis analytical GS energy $E_0$, {dim}")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 10)
plt.savefig(save_dir / "E_N.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
for i, dim in enumerate(dims):
    df = data[('Seq', 'MetroHast', dim)]
    c = colors[i+5]
    ax.plot(df['N'], df['Energy'], color=c, marker='o',
            label=f"MH variational energy $E_T$, {dim}")
    ax.plot(df['N'], E_analyt[dim], 'D', mfc='white', mec=c, ms=6,
            label=f"MH analytical GS $E_0$, {dim}")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 10)
plt.savefig(save_dir / "E_N_MH.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
for i, dim in enumerate(dims):
    df  = data[('Seq', 'MetroHast',    dim)]
    dfN = data[('Seq', 'NumMetroHast', dim)]
    c   = colors[i]

    dim_num = int(dim[0])

    y = df['Energy'] - dfN['Energy'] + dim_num

    ax.plot(df['N'], y, color=c, marker='o',
            label=f"Metropolis, {dim}")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Energy difference $\Delta E + d$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "deltaE.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
for i, dim in enumerate(dims):
    df  = data[('Seq', 'Metro',    dim)]
    dfN = data[('Seq', 'NumMetro', dim)]
    c   = colors[i+3]

    dim_num = int(dim[0])

    y = df['Energy'] - dfN['Energy'] + dim_num

    ax.plot(df['N'], y, color=c, marker='o',
            label=f"MH, {dim}")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Energy difference $\Delta E + d$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "deltaEMH.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))
df1  = data[('Seq', 'Metro', '3D')]
df2 = data[('Seq', 'NumMetro', '3D')]
c1   = colors[2]
c2   = colors[3]

ax.plot(df['N'], df1['Time'], color=c1, marker='o',
        label=f"Metropolis, analytical derivatives, 3D")
ax.plot(df['N'], df2['Time'], color=c2, marker='o',
        label=f"Metropolis, numerical derivatives, 3D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "Time.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))
df1  = data[('Seq', 'MetroHast', '3D')]
df2 = data[('Seq', 'NumMetroHast', '3D')]
c1   = colors[1]
c2   = colors[4]

ax.plot(df['N'], df1['Time'], color=c1, marker='o',
        label=f"MH, analytical derivatives, 3D")
ax.plot(df['N'], df2['Time'], color=c2, marker='o',
        label=f"MH, numerical derivatives, 3D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "TimeMH.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
for i, dim in enumerate(dims):
    df  = data[('Seq', 'Metro',    dim)]
    c   = colors[i]

    ax.plot(df['N'], df['AccRatio'], color=c, marker='o',
            label=f"Metroplois, {dim}")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Acceptance ratio $A$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "AccRate.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))
for i, dim in enumerate(dims):
    df2 = data[('Seq', 'MetroHast', dim)]
    c   = colors[i+3]

    ax.plot(df['N'], df2['AccRatio'], color=c, marker='o',
            label=f"MH, {dim}")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Acceptance ratio $A$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "AccRateMH.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
df  = data[('Seq', 'Metro', '3D')]
dfMPI  = data[('MPI', 'Metro', '3D')]
dfOMP  = data[('OMP', 'Metro', '3D')]

c   = colors[1]
c   = colors[2]

y1 = df['Energy'] - dfOMP['Energy'] + 1
y2 = df['Energy'] - dfMPI['Energy'] + 2

ax.plot(df['N'], y1, color=c1, marker='o',
            label=f"OMP, 3D")

ax.plot(df['N'], y2, color=c2, marker='o',
            label=f"MPI, 3D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Energy difference $\Delta E + b$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "deltaEParallel.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
df  = data[('Seq', 'MetroHast', '3D')]
dfMPI  = data[('MPI', 'MetroHast', '3D')]
dfOMP  = data[('OMP', 'MetroHast', '3D')]

c1   = colors[3]
c2   = colors[4]

y1 = df['Energy'] - dfOMP['Energy'] + 1
y2 = df['Energy'] - dfMPI['Energy'] + 2

ax.plot(df['N'], y1, color=c1, marker='o',
            label=f"OMP, 3D")

ax.plot(df['N'], y2, color=c2, marker='o',
            label=f"MPI, 3D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Energy difference $\Delta E + b$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "deltaEMHParallel.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
df  = data[('Seq', 'Metro', '3D')]
dfMPI  = data[('MPI', 'Metro', '3D')]
dfOMP  = data[('OMP', 'Metro', '3D')]

c1   = colors[1]
c2   = colors[2]
c3   = colors[3]

ax.plot(df['N'], df['Time'], color=c1, marker='o',
            label=f"Serial, 3D")

ax.plot(df['N'], dfOMP['Time'], color=c3, marker='o',
            label=f"OMP, 3D")

ax.plot(df['N'], dfMPI['Time'], color=c2, marker='o',
            label=f"MPI, 3D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "TimeParallel.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))
df  = data[('Seq', 'NumMetro', '3D')]
dfMPI  = data[('MPI', 'NumMetro', '3D')]
dfOMP  = data[('OMP', 'NumMetro', '3D')]

c1   = colors[4]
c2   = colors[5]
c3   = colors[6]

ax.plot(df['N'], df['Time'], color=c1, marker='o',
            label=f"Serial, 3D")

ax.plot(df['N'], dfOMP['Time'], color=c3, marker='o',
            label=f"OMP, 3D")

ax.plot(df['N'], dfMPI['Time'], color=c2, marker='o',
            label=f"MPI, 3D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "TimeParallelNum.pdf")




fig, ax = plt.subplots(figsize=(6.5, 4.6))
df  = data[('Seq', 'MetroHast', '3D')]
dfMPI  = data[('MPI', 'MetroHast', '3D')]
dfOMP  = data[('OMP', 'MetroHast', '3D')]

c1   = colors[1]
c2   = colors[2]
c3   = colors[3]

ax.plot(df['N'], df['Time'], color=c1, marker='o',
            label=f"Serial, 3D")

ax.plot(df['N'], dfOMP['Time'], color=c3, marker='o',
            label=f"OMP, 3D")

ax.plot(df['N'], dfMPI['Time'], color=c2, marker='o',
            label=f"MPI, 3D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "TimeParallelMH.pdf")




fig, ax = plt.subplots(figsize=(6.5, 4.6))
df  = data[('Seq', 'NumMetroHast', '3D')]
dfMPI  = data[('MPI', 'NumMetroHast', '3D')]
dfOMP  = data[('OMP', 'NumMetroHast', '3D')]

c1   = colors[4]
c2   = colors[5]
c3   = colors[6]

ax.plot(df['N'], df['Time'], color=c1, marker='o',
            label=f"Serial, 3D")

ax.plot(df['N'], dfOMP['Time'], color=c3, marker='o',
            label=f"OMP, 3D")

ax.plot(df['N'], dfMPI['Time'], color=c2, marker='o',
            label=f"MPI, 3D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Time $t$, $[s]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "TimeParallelNumMH.pdf")