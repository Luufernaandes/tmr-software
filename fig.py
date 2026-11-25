import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.family':'serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
  'axes.edgecolor':'#888','axes.labelcolor':'#222','xtick.color':'#444','ytick.color':'#444'})
B,O='#257de0','#f47230'
ci=lambda x: stats.t.ppf(.975,len(x)-1)*x.std()/np.sqrt(len(x))
fig,axs=plt.subplots(1,2,figsize=(6.3,2.5),sharey=True)
lab={'base':'Base','dmr':'DMR','tmr':'TMR\nseq.','tmrp':'TMR\nparal.'}
for ax,(f,t) in zip(axs,[('mm','(a) MM, N = 384 (CPU)'),('vs','(b) VS, n = 2$^{24}$ (memória)')]):
    d=pd.read_csv(f'tempo_{f}.csv'); b=d[d.esquema=='base'].t_total.mean()
    xs=np.arange(4)
    for i,s in enumerate(['base','dmr','tmr','tmrp']):
        x=d[d.esquema==s]; r=x.t_replicas.mean()/b; v=x.t_voto.mean()/b; tot=x.t_total.mean()/b
        ax.bar(i,r,0.6,color=B,edgecolor='white',linewidth=1,label='Réplicas (cópia + cálculo)' if i==0 else None)
        if v>0.005: ax.bar(i,v,0.6,bottom=r,color=O,edgecolor='white',linewidth=1,label='Voto / comparação' if i==1 else None)
        ax.errorbar(i,tot,yerr=ci(x.t_total)/b,color='#333',capsize=3,lw=0.8)
        ax.text(i,tot+ci(x.t_total)/b+0.08,f'{tot:.2f}×',ha='center',fontsize=8,color='#222')
    ax.set_xticks(xs,[lab[s] for s in ['base','dmr','tmr','tmrp']]); ax.set_title(t,fontsize=9)
    ax.yaxis.grid(True,color='#e5e5e5',lw=0.6); ax.set_axisbelow(True); ax.set_ylim(0,4.8)
axs[0].set_ylabel('Tempo normalizado (Base = 1)')
h,l=axs[1].get_legend_handles_labels(); h0,l0=axs[0].get_legend_handles_labels()
fig.legend(h,l,loc='upper center',ncol=2,frameon=False,fontsize=8,bbox_to_anchor=(0.5,1.04))
fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig('fig_tempo.pdf'); fig.savefig('fig_tempo.png',dpi=150)

v=pd.read_csv('votador.csv')
fig,ax=plt.subplots(figsize=(3.3,2.2))
for j,(vt,c,n) in enumerate([('palavra',B,'Por palavra (com if/else)'),('bit_a_bit',O,'Bit a bit (sem desvios condicionais)')]):
    for i,sc in enumerate(['limpo','ruidoso']):
        x=v[(v.cenario==sc)&(v.votador==vt)].ns_por_palavra; m=x.mean()
        ax.bar(i+(j-0.5)*0.36,m,0.34,color=c,edgecolor='white',label=n if i==0 else None)
        ax.errorbar(i+(j-0.5)*0.36,m,yerr=ci(x),color='#333',capsize=2,lw=0.8)
        ax.text(i+(j-0.5)*0.36,m+ci(x)+0.25,f'{m:.1f}',ha='center',fontsize=8)
ax.set_xticks([0,1],['Réplicas iguais','50% das palavras\ncom 1 divergência']); ax.set_ylabel('ns por palavra')
ax.yaxis.grid(True,color='#e5e5e5',lw=0.6); ax.set_axisbelow(True); ax.set_ylim(0,13.5)
ax.legend(frameon=False,fontsize=7.5,loc='upper left')
fig.tight_layout(); fig.savefig('fig_votador.pdf'); fig.savefig('fig_votador.png',dpi=150)
