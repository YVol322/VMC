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

save_dir = Path('Results/Figures/Fermions/Interacting/Program1')
save_dir.mkdir(parents=True, exist_ok=True)


base = Path('Results/Tables/Fermions/Interacting/Program1')
files = {
    "J2"     : base / 'OMP_eta_h_N=2.cvv',
    "J2MH"  : base / 'OMP_eta_h_N=2_MH.cvv',
    "PJ2"     : base / 'OMP_eta_h_N=2_PJ.cvv',
    "PJ2MH"  : base / 'OMP_eta_h_N=2_MH_PJ.cvv',
}

cols = ['Iter', 'Eta', 'H', 'Energy', 'grad', 'beta']

data = {}
for key, path in files.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path, sep='\t')
    data[key] = df[cols]



df = data['J2'] 

pivot = df.pivot(index='Eta', columns='H', values='Iter')
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
plt.savefig(save_dir / 'Ferm_J2_heatmap_iter.pdf')


pivot = df.pivot(index='Eta', columns='H', values='Energy')

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
plt.savefig(save_dir / 'Ferm_J2_heatmap_En.pdf')


pivot = df.pivot(index='Eta', columns='H', values='beta')

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational parameter $\beta$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Variational parameter $\beta$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_J2_heatmap_beta_MH.pdf')




df = data['J2MH'] 

pivot = df.pivot(index='Eta', columns='H', values='Iter')
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
plt.savefig(save_dir / 'Ferm_J2_heatmap_iter_MH.pdf')



pivot = df.pivot(index='Eta', columns='H', values='Energy')

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
plt.savefig(save_dir / 'Ferm_J2_heatmap_En_MH.pdf')



pivot = df.pivot(index='Eta', columns='H', values='beta')

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational parameter $\beta$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Variational parameter $\beta$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_J2_heatmap_beta_MH.pdf')


df = data['PJ2'] 

pivot = df.pivot(index='Eta', columns='H', values='Iter')
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
plt.savefig(save_dir / 'Ferm_PJ2_heatmap_iter.pdf')


pivot = df.pivot(index='Eta', columns='H', values='Energy')

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
plt.savefig(save_dir / 'Ferm_PJ2_heatmap_En.pdf')



pivot = df.pivot(index='Eta', columns='H', values='beta')

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational parameter $\beta$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Variational parameter $\beta$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_PJ2_heatmap_beta.pdf')


df = data['PJ2MH'] 

pivot = df.pivot(index='Eta', columns='H', values='Iter')
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
plt.savefig(save_dir / 'Ferm_PJ2_heatmap_iter_MH.pdf')


pivot = df.pivot(index='Eta', columns='H', values='Energy')

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
plt.savefig(save_dir / 'Ferm_PJ2_heatmap_En_MH.pdf')



pivot = df.pivot(index='Eta', columns='H', values='beta')

plt.figure(figsize=(6.5, 4.6))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational parameter $\beta$'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'Variational parameter $\beta$', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_PJ2_heatmap_beta_MH.pdf')


c = sns.color_palette("colorblind")
plt.style.use("ggplot")

base4 = Path('Results/Tables/Fermions/Interacting/Program2')
fname = 'OMP_J_MH.cvv'
path  = base4 / fname

cols = ['Iter', 'Energy', 'Grad', 'Beta']

if not path.exists():
    raise FileNotFoundError(f"Missing file: {path}")
df = pd.read_csv(path, sep='\t')[cols]



fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(df['Iter'], df['Energy'], color=c[0], marker='o',
        label=f"Jastow MH")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r'Variational energy $E_T$ $[a.u.]$', fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "J2_En.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.0))

ax.plot(df['Iter'], df['Grad'], color=c[1], marker='o',
        label=f"Jastow MH")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\beta_{12}}\Psi_T\rangle$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "J2_Grad.pdf")
