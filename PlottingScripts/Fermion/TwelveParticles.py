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

save_dir = Path('Results/Figures/Fermions/Interacting/Program4')
save_png = Path('Results/Figures/PNG/Fermions/Interacting/Program4')
save_dir.mkdir(parents=True, exist_ok=True)
save_png.mkdir(parents=True, exist_ok=True)


base = Path('Results/Tables/Fermions/Interacting/Program5')
files = {
    "J12"     : base / 'OMP_eta_h_N=12.cvv',
    "J12MH"  : base / 'OMP_eta_h_N=12_MH.cvv',
    "PJ12"     : base / 'OMP_eta_h_N=12_PJ.cvv',
    "PJ12MH"  : base / 'OMP_eta_h_N=12_MH_PJ.cvv',
}

cols = ['Iter', 'Eta', 'H', 'Energy', 'grad']

data = {}
for key, path in files.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path, sep='\t')
    data[key] = df[cols]



df = data['J12'] 

pivot = df.pivot(index='Eta', columns='H', values='Iter')
pivot = pivot.replace([np.inf, -np.inf], 0)
pivot = pivot.fillna(0)


plt.figure(figsize=(6.5, 5.1))
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
plt.savefig(save_dir / 'Ferm_J12_heatmap_iter.pdf')
plt.savefig(save_png / 'Ferm_J12_heatmap_iter.png')


pivot = df.pivot(index='Eta', columns='H', values='Energy')
pivot = pivot.replace([np.inf, -np.inf], 0)
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 5.1))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational Energy $E_T$'}
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
plt.savefig(save_dir / 'Ferm_J12_heatmap_En.pdf')
plt.savefig(save_png / 'Ferm_J12_heatmap_En.png')


pivot = df.pivot(index='Eta', columns='H', values='grad')
pivot = pivot.replace([np.inf, -np.inf], 100)
pivot = pivot.fillna(100)


plt.figure(figsize=(6.5, 5.1))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".1e",
    cbar_kws={'label': r'L2 norm'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'L2 norm', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_J12_heatmap_grad_MH.pdf')
plt.savefig(save_png / 'Ferm_J12_heatmap_grad_MH.png')




df = data['J12MH'] 

pivot = df.pivot(index='Eta', columns='H', values='Iter')
pivot = pivot.replace([np.inf, -np.inf], 0)
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 5.1))
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
plt.savefig(save_dir / 'Ferm_J12_heatmap_iter_MH.pdf')
plt.savefig(save_png / 'Ferm_J12_heatmap_iter_MH.png')



pivot = df.pivot(index='Eta', columns='H', values='Energy')
pivot = pivot.replace([np.inf, -np.inf], 0)
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 5.1))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational Energy $E_T$'}
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
plt.savefig(save_dir / 'Ferm_J12_heatmap_En_MH.pdf')
plt.savefig(save_png / 'Ferm_J12_heatmap_En_MH.png')



pivot = df.pivot(index='Eta', columns='H', values='grad')
pivot = pivot.replace([np.inf, -np.inf], 100)
pivot = pivot.fillna(100)

plt.figure(figsize=(6.5, 5.1))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".1e",
    cbar_kws={'label': r'L2 norm'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'L2 norm', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_J12_heatmap_grad_MH.pdf')
plt.savefig(save_png / 'Ferm_J12_heatmap_grad_MH.png')


df = data['PJ12'] 

pivot = df.pivot(index='Eta', columns='H', values='Iter')
pivot = pivot.replace([np.inf, -np.inf], 0)
pivot = pivot.fillna(0)
plt.figure(figsize=(6.5, 5.1))
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
plt.savefig(save_dir / 'Ferm_PJ12_heatmap_iter.pdf')
plt.savefig(save_png / 'Ferm_PJ12_heatmap_iter.png')


pivot = df.pivot(index='Eta', columns='H', values='Energy')
pivot = pivot.replace([np.inf, -np.inf], 0)
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 5.1))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational Energy $E_T$'}
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
plt.savefig(save_dir / 'Ferm_PJ12_heatmap_En.pdf')
plt.savefig(save_png / 'Ferm_PJ12_heatmap_En.png')



pivot = df.pivot(index='Eta', columns='H', values='grad')
pivot = pivot.replace([np.inf, -np.inf], 100)
pivot = pivot.fillna(100)

plt.figure(figsize=(6.5, 5.1))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".1e",
    cbar_kws={'label': r'L2 norm'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'L2 norm', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_PJ12_heatmap_grad.pdf')
plt.savefig(save_png / 'Ferm_PJ12_heatmap_grad.png')


df = data['PJ12MH'] 

pivot = df.pivot(index='Eta', columns='H', values='Iter')
pivot = pivot.replace([np.inf, -np.inf], 0)
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 5.1))
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
plt.savefig(save_dir / 'Ferm_PJ12_heatmap_iter_MH.pdf')
plt.savefig(save_png / 'Ferm_PJ12_heatmap_iter_MH.png')


pivot = df.pivot(index='Eta', columns='H', values='Energy')
pivot = pivot.replace([np.inf, -np.inf], 0)
pivot = pivot.fillna(0)

plt.figure(figsize=(6.5, 5.1))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".3f",
    cbar_kws={'label': r'Variational Energy $E_T$'}
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
plt.savefig(save_dir / 'Ferm_PJ12_heatmap_En_MH.pdf')
plt.savefig(save_png / 'Ferm_PJ12_heatmap_En_MH.png')



pivot = df.pivot(index='Eta', columns='H', values='grad')
pivot = pivot.replace([np.inf, -np.inf], 100)
pivot = pivot.fillna(100)

plt.figure(figsize=(6.5, 5.1))
ax = sns.heatmap(
    pivot,
    cmap='viridis',
    annot=True,
    fmt=".1e",
    cbar_kws={'label': r'L2 norm'}
)
ax.set_xlabel(r'Step size $h$', fontsize=16)
ax.set_ylabel(r'Learning rate $\eta$', fontsize=16)
# Set tick label sizes
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
# Adjust colorbar label and tick sizes
cbar = ax.collections[0].colorbar
cbar.set_label(r'L2 norm', fontsize=16)
cbar.ax.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig(save_dir / 'Ferm_PJ12_heatmap_grad_MH.pdf')
plt.savefig(save_png / 'Ferm_PJ12_heatmap_grad_MH.png')
#plt.show()

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
ax.set_ylabel(r"Variational energy $E_T$", fontsize = 16)
plt.legend(fontsize = 12)
plt.savefig(save_dir / "J12_En.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(df['Iter'], df['Grad'], color=c[1], marker='o',
        label=f"Jastow MH")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\beta_{12}}\Psi_T\rangle$", fontsize = 16)
plt.legend(fontsize = 12)
plt.savefig(save_dir / "J12_Grad.pdf")
