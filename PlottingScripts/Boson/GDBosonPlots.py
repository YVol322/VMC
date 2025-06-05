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

save_dir = Path('Results/Figures/Bosons/Program2')
save_dir.mkdir(parents=True, exist_ok=True)

colors = sns.color_palette("colorblind")
plt.style.use("ggplot")


base4 = Path('Results/Tables/Bosons/Program4')


files4 = {
    'SeqGrad2D':    'N_Grad_El_2D.cvv',
    'MH_SeqGrad2D': 'MH_N_Grad_El_2D.cvv',
}


cols = ['N', 'Grad', 'Energy']


data4 = {}
for key, fname in files4.items():
    path = base4 / fname
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path, sep='\t')
    data4[key] = df[cols]




fig, ax = plt.subplots(figsize=(6.5, 4.6))
df = data4['SeqGrad2D']
dfMH = data4['MH_SeqGrad2D']
c1 = colors[0]
c2 = colors[1]

ax.plot(df['N'], df['Energy'], color=c1, marker='o',
        label=f"Metropolis, 2D")
ax.plot(df['N'], dfMH['Energy'],'--', color=c2,
        label=f"MH, 2D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Variational energy $E_T$ $[a.u.]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "GD_E_N.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))
df = data4['SeqGrad2D']
dfMH = data4['MH_SeqGrad2D']
c1 = colors[2]
c2 = colors[3]

ax.plot(df['N'], df['Grad'], color=c1, marker='o',
        label=f"Metropolis, 2D")
ax.plot(df['N'], dfMH['Grad'], '--', color=c2,
        label=f"MH, 2D")

ax.set_xlabel(r"Number of particles $N$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "GD_Grad_N.pdf")



base5 = Path('Results/Tables/Bosons/Program5')

files5 = {
    'N2_MH':   '2D_N=2_alpha=0.1_eta=1e-2_MH.cvv',
    'N2':      '2D_N=2_alpha=0.1_eta=1e-2.cvv',
    'N10_MH':  '2D_N=10_alpha=0.1_eta=1e-2_MH.cvv',
    'N10':     '2D_N=10_alpha=0.1_eta=1e-2.cvv',
}


cols = ['Iter', 'Energy', 'Grad', 'Alpha']

data5 = {}
for key, fname in files5.items():
    path = base5 / fname
    if not path.exists():
        raise FileNotFoundError(f"Could not find {path}")
    df = pd.read_csv(path, sep='\t')
    
    data5[key] = df[cols]


df = data5['N2']
dfMH = data5['N2_MH']
df2   = df[df['Iter'] > 2]
dfMH2 = dfMH[dfMH['Iter'] > 2]
c1 = colors[0]
c2 = colors[1]
c3 = colors[2]
c4 = colors[3]
c5 = colors[4]
c6 = colors[5]
c7 = colors[6]

fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(df2['Iter'], df2['Grad'], color=c1,
        label=f"Metropolis, $N=2$, 2D")
ax.plot(dfMH2['Iter'], dfMH2['Grad'], '--', color=c2,
        label=f"MH, $N=2$, 2D")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "GD_Grad_i2.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=2.0, color=c1, linewidth=1.5, label="Analytical GS energy $E_0=2$", linestyle=':')

ax.plot(df2['Iter'], df2['Energy'],  color=c1,
        label=r"Metropolis variational energy $E_T$, $N=2$, 2D")
ax.plot(dfMH2['Iter'], dfMH2['Energy'], '--', color=c2,
        label=r"MH variational energy $E_T$, $N=2$, 2D")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ $[a.u.]$", fontsize = 16)
ax.legend(fontsize = 10)
plt.savefig(save_dir / "GD_En_i2.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=0.5, color=c1, linewidth=1.5, label=r"Optimal $\alpha=0.5$", linestyle = ':')
ax.plot(df2['Iter'], df2['Alpha'],  color=c1,
        label=r"Metropolis variational parameter $\alpha$, $N=2$, 2D")
ax.plot(dfMH2['Iter'], dfMH2['Alpha'], '--', color=c2,
        label=r"MH variational parameter $\alpha$, $N=2$, 2D")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Variational parameter $\alpha$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "GD_alp_i2.pdf")



df = data5['N10']
dfMH = data5['N10_MH']
df2   = df[df['Iter'] > 2]
dfMH2 = dfMH[dfMH['Iter'] > 2]

ig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.plot(df2['Iter'], df2['Grad'], color=c3,
        label=f"Metropolis, $N=10$, 2D")
ax.plot(dfMH2['Iter'], dfMH2['Grad'], '--', color=c4,
        label=f"MH, $N=10$, 2D")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "GD_Grad_i10.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=10.0, color=c1, linewidth=1.5, label="Analytical GS energy $E_0=10$", linestyle = ":")
ax.plot(df2['Iter'], df2['Energy'],  color=c3,
        label=r"Metropolis variational energy $E_T$, $N=10$, 2D")
ax.plot(dfMH2['Iter'], dfMH2['Energy'], '--', color=c4,
        label=r"MH variational energy $E_T$, $N=10$, 2D")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Energy $E$ $[a.u.]$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "GD_En_i10.pdf")


fig, ax = plt.subplots(figsize=(6.5, 4.6))

ax.axhline(y=0.5, color=c1, linewidth=1.5, label=r"Optimal $\alpha=0.5$", linestyle = ":")
ax.plot(df2['Iter'], df2['Alpha'],  color=c3,
        label=r"Metropolis variational parameter $\alpha$, $N=10$, 2D")
ax.plot(dfMH2['Iter'], dfMH2['Alpha'], '--', color=c4,
        label=r"MH variational parameter $\alpha$, $N=2$, 2D")

ax.set_xlabel(r"Number of iterations $i$", fontsize = 16)
ax.set_ylabel(r"Variational parameter $\alpha$", fontsize = 16)
ax.legend(fontsize = 12)
plt.savefig(save_dir / "GD_alp_i10.pdf")




base = Path('Results/Tables/Bosons/Program6')
files = {
    "Metropolis"     : base / '2D_N=5_alpha=1eta_h.cvv',
    "Metropolis-Hastings"  : base / '2D_N=5_alpha=1eta_h_MH.cvv'
}

df_list = []
for method, path in files.items():
    df = pd.read_csv(path, sep='\t')
    df['Method'] = method
    df_list.append(df)

df = pd.concat(df_list, ignore_index=True)


for method in df['Method'].unique():
    sub = df[df['Method'] == method]
    pivot = sub.pivot_table(
        index='Eta', columns='H', values='Iter', aggfunc='mean'
    )
    
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
    plt.xlabel(r'Step size $h$', fontsize = 16)
    plt.ylabel(r'Learning rate $\eta$', fontsize = 16)
	
    plt.savefig(save_dir / f'heatmap_{method}.pdf')
    