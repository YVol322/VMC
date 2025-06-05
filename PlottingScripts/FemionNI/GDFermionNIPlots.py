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

save_dir = Path('Results/Figures/Fermions/NonInteracting/Program2')
save_dir.mkdir(parents=True, exist_ok=True)


base4 = Path('Results/Tables/Fermions/NonInteracting/Program4')
fname = 'OMP_GD.cvv'
path  = base4 / fname

cols = ['Energy', 'Time', 'Grad']

if not path.exists():
    raise FileNotFoundError(f"Missing file: {path}")
df = pd.read_csv(path, sep='\t')[cols]


c = sns.color_palette("colorblind")
plt.style.use("ggplot")

nparts = np.array([2, 6, 12])

EnergyMetro2 = df['Energy'][0]
EnergyMetro6 = df['Energy'][2]
EnergyMetro12 = df['Energy'][4]

EnergyMH2 = df['Energy'][1]
EnergyMH6 = df['Energy'][3]
EnergyMH12 = df['Energy'][5]

EnMetroarr = np.array([EnergyMetro2, EnergyMetro6, EnergyMetro12])
EnMHarr = np.array([EnergyMH2, EnergyMH6, EnergyMH12])

GradMetro2 = df['Grad'][0]
GradMetro6 = df['Grad'][2]
GradMetro12 = df['Grad'][4]

GradMH2 = df['Grad'][1]
GradMH6 = df['Grad'][3]
GradMH12 = df['Grad'][5]

GradMetroArr = np.array([GradMetro2, GradMetro6, GradMetro12])
GradMHArr = np.array([GradMH2, GradMH6, GradMH12])


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, EnMetroarr, color=c[0], marker='o',
        label=f"Metropolis")
ax.plot(nparts, EnMHarr,'--', color=c[1], marker='o',
        label=f"MH")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Variational energy $E_T$ [a.u.]", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "FermGD_E_N.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(nparts, GradMetroArr, color=c[2], marker='o',
        label=f"Metropolis")
ax.plot(nparts, GradMHArr, '--', color=c[3], marker='o',
        label=f"MH")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
#
plt.savefig(save_dir / "FermGD_Grad_N.pdf")


# Directory and filenames for Program5
base5 = Path('Results/Tables/Fermions/NonInteracting/Program5')

files5 = {
    'N2_MH':   'N=2_alpha=0.1_h=0.01_eta=0.1_MH.cvv',
    'N2':      'N=2_alpha=0.1_h=1_eta=0.1.cvv',
    'N6_MH':  'N=6_alpha=0.1_h=0.01_eta=0.1_MH.cvv',
    'N6':     'N=6_alpha=0.1_h=0.1_eta=0.1.cvv',
    'N12_MH':  'N=12_alpha=0.1_h=0.001_eta=0.1_MH.cvv',
    'N12':     'N=12_alpha=0.1_h=0.1_eta=0.1.cvv',
}

# assume these columns exist in all four
cols = ['Iter', 'Energy', 'Grad', 'Alpha']

data5 = {}
for key, fname in files5.items():
    path = base5 / fname
    if not path.exists():
        raise FileNotFoundError(f"Could not find {path}")
    df = pd.read_csv(path, sep='\t')
    # if you only need a subset of columns:
    data5[key] = df[cols]



fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(data5['N2']['Iter'], data5['N2']['Grad'], color=c[0],
        label=f"Metropolis $N=2$")
ax.plot(data5['N2_MH']['Iter'], data5['N2_MH']['Grad'], '--', color=c[1],
        label=f"MH $N=2$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_Grad_i2.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(data5['N6']['Iter'], data5['N6']['Grad'], color=c[0],
        label=f"Metropolis $N=6$")
ax.plot(data5['N6_MH']['Iter'], data5['N6_MH']['Grad'], '--', color=c[1],
        label=f"MH $N=6$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_Grad_i6.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(data5['N12']['Iter'], data5['N12']['Grad'], color=c[2],
        label=f"Metropolis $N=12$")
ax.plot(data5['N12_MH']['Iter'], data5['N12_MH']['Grad'], '--', color=c[3],
        label=f"MH $N=12$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_Grad_i12.pdf")



fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=2.0, color=c[0], linewidth=1.5, label="Analytical GS energy $E_0=2$", linestyle = ':')

ax.plot(data5['N2']['Iter'], data5['N2']['Energy'],  color=c[1],
        label=r"Metropolis variational energy $E_T$, $N=2$")
ax.plot(data5['N2_MH']['Iter'], data5['N2_MH']['Energy'], '--', color=c[2],
        label=r"MH variational energy $E_T$, $N=2$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ $[a.u.]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_En_i2.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=10.0, color=c[0], linewidth=1.5, label="Analytical GS energy $E_0=10$", linestyle = ':')

ax.plot(data5['N6']['Iter'], data5['N6']['Energy'],  color=c[3],
        label=r"Metropolis variational energy $E_T$, $N=6$")
ax.plot(data5['N6_MH']['Iter'], data5['N6_MH']['Energy'], '--', color=c[4],
        label=r"MH variational energy $E_T$, $N=6$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ $[a.u.]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_En_i6.pdf")

fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=28.0, color=c[0], linewidth=1.5, label="Analytical GS energy $E_0=28$", linestyle = ':')

ax.plot(data5['N12']['Iter'], data5['N12']['Energy'],  color=c[5],
        label=r"Metropolis variational energy $E_T$, $N=12$")
ax.plot(data5['N12_MH']['Iter'], data5['N12_MH']['Energy'], '--', color=c[6],
        label=r"Metropolis variational energy $E_T$, $N=12$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ $[a.u.]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_En_i12.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=0.5, color=c[0], linewidth=1.5, label=r"Optimal $\alpha=0.5$", linestyle = ":")
ax.plot(data5['N2']['Iter'], data5['N2']['Alpha'],  color=c[1],
        label=r"Metropolis variational parameter $\alpha$, $N=2$")
ax.plot(data5['N2_MH']['Iter'], data5['N2_MH']['Alpha'], '--', color=c[2],
        label=r"MH variational parameter $\alpha$, $N=2$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Variational parameter $\alpha$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_alp_i2.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=0.5, color=c[0], linewidth=1.5, label=r"Optimal $\alpha=0.5$", linestyle = ":")
ax.plot(data5['N6']['Iter'], data5['N6']['Alpha'],  color=c[3],
        label=f"Metropolis $N=6$")
ax.plot(data5['N6_MH']['Iter'], data5['N6_MH']['Alpha'], '--', color=c[4],
        label=f"MH $N=6$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Variational parameter $\alpha$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_alp_6.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=0.5, color=c[0], linewidth=1.5, label=r"Optimal $\alpha=0.5$", linestyle = ":")
ax.plot(data5['N12']['Iter'], data5['N12']['Alpha'],  color=c[5],
        label=f"Metropolis $N=12$")
ax.plot(data5['N12_MH']['Iter'], data5['N12_MH']['Alpha'], '--', color=c[6],
        label=f"MH $N=12$")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Variational parameter $\alpha$", fontsize = 16)
ax.legend(fontsize = 12)
plt.tight_layout()
plt.savefig(save_dir / "Ferm_GD_alp_i12.pdf")



base = Path('Results/Tables/Fermions/NonInteracting/Program6')
files = {
    "Me2"     : base / 'OMP_eta_h_N=2.cvv',
    "MH2"  : base / 'OMP_eta_h_N=2_MH.cvv',
    "Me6"     : base / 'OMP_eta_h_N=6.cvv',
    "MH6"  : base / 'OMP_eta_h_N=6_MH.cvv',
    "Me12"     : base / 'OMP_eta_h_N=12.cvv',
    "MH12"  : base / 'OMP_eta_h_N=12_MH.cvv',
}

cols = ['Iter', 'Eta', 'H', 'Energy']

data = {}
for key, path in files.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path, sep='\t')
    data[key] = df[cols]



df = data['Me12']  # or whatever key

# Pivot: rows=Eta, cols=H, values=Iter
pivot = df.pivot(index='Eta', columns='H', values='Iter')
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".0f",
    cbar_kws={'label': r'Number of iterations $i$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Number of iterations $i$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_heatmap_iter12.pdf')


df = data['MH12']  # or whatever key
# Pivot: rows=Eta, cols=H, values=Iter
pivot = df.pivot(index='Eta', columns='H', values='Iter')
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".0f",
    cbar_kws={'label': r'Number of iterations $i$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Number of iterations $i$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()

#
plt.savefig(save_dir / 'Ferm_heatmap_iter12_MH.pdf')


df = data['Me12']  # or whatever key

# Pivot: rows=Eta, cols=H, values=Iter
pivot = df.pivot(index='Eta', columns='H', values='Energy')
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational energy $E_T$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Variational energy $E_T$ $[a.u.]$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()

#
plt.savefig(save_dir / 'Ferm_heatmap_En12.pdf')



df = data['MH12']  # or whatever key
# Pivot: rows=Eta, cols=H, values=Iter
pivot = df.pivot(index='Eta', columns='H', values='Energy')
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational energy $E_T$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Variational energy $E_T$ $[a.u.]$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_heatmap_En12_MH.pdf')



df = data['Me6']  # or whatever key

# Pivot: rows=Eta, cols=H, values=Iter
pivot = df.pivot(index='Eta', columns='H', values='Iter')
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".0f",
    cbar_kws={'label': r'Number of iterations $i$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Number of iterations $i$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_heatmap_iter6.pdf')


df = data['MH6']  # or whatever key
# Pivot: rows=Eta, cols=H, values=Iter
pivot = df.pivot(index='Eta', columns='H', values='Iter')
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".0f",
    cbar_kws={'label': r'Number of iterations $i$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Number of iterations $i$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_heatmap_iter6_MH.pdf')


df = data['Me6']  # or whatever key

# Pivot: rows=Eta, cols=H, values=Iter
pivot = df.pivot(index='Eta', columns='H', values='Energy')
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational energy $E_T$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Variational energy $E_T$ $[a.u.]$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_heatmap_En6.pdf')


df = data['MH6']  # or whatever key
# Pivot: rows=Eta, cols=H, values=Iter
pivot = df.pivot(index='Eta', columns='H', values='Energy')
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational energy $E_T$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Variational energy $E_T$ $[a.u.]$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_heatmap_En6_MH.pdf')
