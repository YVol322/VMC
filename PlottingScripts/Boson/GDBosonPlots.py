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
#
## Directory containing your new tables
#base4 = Path('Results/Tables/Bosons/Program4')
#
## The two files you want to load, with short keys
#files4 = {
#    'SeqGrad2D':    'N_Grad_El_2D.cvv',
#    'MH_SeqGrad2D': 'MH_N_Grad_El_2D.cvv',
#}
#
## Columns we care about (adjust if the files have different names)
#cols = ['N', 'Grad', 'Energy']
#
## Read them into a dict of DataFrames
#data4 = {}
#for key, fname in files4.items():
#    path = base4 / fname
#    if not path.exists():
#        raise FileNotFoundError(f"Missing file: {path}")
#    df = pd.read_csv(path, sep='\t')
#    data4[key] = df[cols]
#
#
#colors = sns.color_palette("colorblind")
#plt.style.use("ggplot")



#fig, ax = plt.subplots(figsize=(6.5, 4.0))
#df = data4['SeqGrad2D']
#dfMH = data4['MH_SeqGrad2D']
#c1 = colors[0]
#c2 = colors[1]
#
#ax.plot(df['N'], df['Energy'], color=c1, marker='o',
#        label=f"2D Metropolis")
#ax.plot(df['N'], dfMH['Energy'],'--', color=c2, marker='o',
#        label=f"2D MH")
#
#ax.set_xlabel(r"Number of particles $N$")
#ax.set_ylabel(r"GS energy $E_0$")
#plt.legend()
##plt.savefig(save_dir / "GD_E_N.pdf")
#
#
#fig, ax = plt.subplots(figsize=(6.5, 4.0))
#df = data4['SeqGrad2D']
#dfMH = data4['MH_SeqGrad2D']
#c1 = colors[0]
#c2 = colors[1]
#
#ax.plot(df['N'], df['Grad'], color=c1, marker='o',
#        label=f"2D Metropolis")
#ax.plot(df['N'], dfMH['Grad'], '--', color=c2, marker='o',
#        label=f"2D MH")
#
#ax.set_xlabel(r"Number of particles $N$")
#ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$")
#plt.legend()
#plt.savefig(save_dir / "GD_Grad_N.pdf")


## Directory and filenames for Program5
#base5 = Path('Results/Tables/Bosons/Program5')
#
#files5 = {
#    'N2_MH':   '2D_N=2_alpha=0.1_eta=1e-2_MH.cvv',
#    'N2':      '2D_N=2_alpha=0.1_eta=1e-2.cvv',
#    'N10_MH':  '2D_N=10_alpha=0.1_eta=1e-2_MH.cvv',
#    'N10':     '2D_N=10_alpha=0.1_eta=1e-2.cvv',
#}
#
## assume these columns exist in all four
#cols = ['Iter', 'Energy', 'Grad', 'Alpha']
#
#data5 = {}
#for key, fname in files5.items():
#    path = base5 / fname
#    if not path.exists():
#        raise FileNotFoundError(f"Could not find {path}")
#    df = pd.read_csv(path, sep='\t')
#    # if you only need a subset of columns:
#    data5[key] = df[cols]
#
#
#df = data5['N2']
#dfMH = data5['N2_MH']
#df2   = df[df['Iter'] > 2]
#dfMH2 = dfMH[dfMH['Iter'] > 2]
#c1 = colors[0]
#c2 = colors[1]
#c3 = colors[2]
#c4 = colors[3]
#c5 = colors[4]
#c6 = colors[5]
#c7 = colors[6]
#
#fig, ax = plt.subplots(figsize=(6.5, 4.0))
#
#ax.plot(df2['Iter'], df2['Grad'], color=c1,
#        label=f"2D Metropolis $N=2$")
#ax.plot(dfMH2['Iter'], dfMH2['Grad'], '--', color=c2,
#        label=f"2D MH $N=2$")
#
#ax.set_xlabel(r"Iteration $i$")
#ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$")
#plt.legend()
#plt.savefig(save_dir / "GD_Grad_i2.pdf")
#
#
#fig, ax = plt.subplots(figsize=(6.5, 4.0))
#
#ax.axhline(y=2.0, color=c1, linewidth=1.5, label="Analytical $E=2$")
#
#ax.plot(df2['Iter'], df2['Energy'], ':',  color=c3,
#        label=f"2D Metropolis $N=2$")
#ax.plot(dfMH2['Iter'], dfMH2['Energy'], '--', color=c4,
#        label=f"2D MH $N=2$")
#
#ax.set_xlabel(r"Iteration $i$")
#ax.set_ylabel(r"GS Energy $E_0$")
#plt.legend()
#plt.savefig(save_dir / "GD_En_i2.pdf")
#
#
#fig, ax = plt.subplots(figsize=(6.5, 4.0))
#
#ax.axhline(y=0.5, color=c1, linewidth=1.5, label=r"Analytical $\alpha=0.5$")
#ax.plot(df2['Iter'], df2['Alpha'], ':',  color=c5,
#        label=f"2D Metropolis $N=2$")
#ax.plot(dfMH2['Iter'], dfMH2['Alpha'], '--', color=c6,
#        label=f"2D MH $N=2$")
#
#ax.set_xlabel(r"Iteration $i$")
#ax.set_ylabel(r"Variational parameter $\alpha$")
#plt.legend()
#plt.savefig(save_dir / "GD_alp_i2.pdf")
#
#
#
#df = data5['N10']
#dfMH = data5['N10_MH']
#df2   = df[df['Iter'] > 2]
#dfMH2 = dfMH[dfMH['Iter'] > 2]
#
#ig, ax = plt.subplots(figsize=(6.5, 4.0))
#
#ax.plot(df2['Iter'], df2['Grad'], color=c1,
#        label=f"2D Metropolis $N=10$")
#ax.plot(dfMH2['Iter'], dfMH2['Grad'], '--', color=c2,
#        label=f"2D MH")
#
#ax.set_xlabel(r"Iteration $i$")
#ax.set_ylabel(r"Gradient $\langle \nabla_{\alpha}\Psi_T\rangle$")
#plt.legend()
#plt.savefig(save_dir / "GD_Grad_i10.pdf")
#
#
#fig, ax = plt.subplots(figsize=(6.5, 4.0))
#
#ax.axhline(y=10.0, color=c1, linewidth=1.5, label="Analytical $E=10$")
#ax.plot(df2['Iter'], df2['Energy'],  color=c3,
#        label=f"2D Metropolis $N=10$")
#ax.plot(dfMH2['Iter'], dfMH2['Energy'], '--', color=c4,
#        label=f"2D MH")
#
#ax.set_xlabel(r"Iteration $i$")
#ax.set_ylabel(r"GS Energy $E_0$")
#plt.legend()
#plt.savefig(save_dir / "GD_En_i10.pdf")
#
#
#fig, ax = plt.subplots(figsize=(6.5, 4.0))
#
#ax.axhline(y=0.5, color=c1, linewidth=1.5, label=r"Analytical $\alpha=0.5$")
#ax.plot(dfMH2['Iter'], dfMH2['Alpha'], color=c2,
#        label=f"2D MH")
#ax.plot(df2['Iter'], df2['Alpha'], ':',  color=c3,
#        label=f"2D Metropolis $N=10$")
#
#ax.set_xlabel(r"Iteration $i$")
#ax.set_ylabel(r"Variational parameter $\alpha$")
#plt.legend()
#plt.savefig(save_dir / "GD_alp_i10.pdf")




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
    
    plt.figure(figsize=(6, 4.4))
    sns.heatmap(
        pivot,
        cmap='viridis',
        annot=True, 
        fmt=".0f", 
        cbar_kws={'label': r'Iteration $i$'}
    )
    plt.xlabel(r'Step size $h$')
    plt.ylabel(r'Learning rate $\eta$')
	# or, if you still need more:
    plt.savefig(save_dir / f'heatmap_{method}.pdf')
    #plt.show()