from pathlib import Path
import os
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(os.environ.get('ANALYSIS_OUTPUT_DIR','outputs')).resolve();A=O/'Analysis';W=Path(os.environ.get('PRIVATE_WORK_DIR','work/private'));W.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','axes.linewidth':.7,'legend.frameon':False})
blue='#24658B';orange='#BC633A';gray='#525B64'
def save(fig,n):
 fig.savefig(O/f'Supplementary_Figure_S{1 if n==1 else 2}.tif',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'});fig.savefig(O/f'Supplementary_Figure_S{1 if n==1 else 2}.svg',bbox_inches='tight');fig.savefig(W/f'Figure_{n}.png',dpi=160,bbox_inches='tight');plt.close(fig)
x=pd.read_csv(A/'primary_monthly_curves.csv');fig,ax=plt.subplots(1,2,figsize=(7.2,3.2));p=x[x.method=='I'];q=x[x.method=='II'];t=p.Interval.to_numpy()/12
ax[0].plot(t,p['Observed Survival (Cum)'],color=gray,label='Observed');ax[0].plot(t,p['Expected Survival (Cum)'],color=blue,label='Expected Ederer I');ax[0].plot(t,q['Expected Survival (Cum)'],color=orange,ls='--',label='Expected Ederer II');ax[0].set(xlabel='Years since diagnosis',ylabel='Survival probability',ylim=(.68,1.01));ax[0].legend(fontsize=8,loc='lower left')
for method,col,ls in [('I',blue,'-'),('II',orange,'--')]:
 z=x[x.method==method];v=z['Expected Survival (Cum)'].to_numpy()-z['Observed Survival (Cum)'].to_numpy();integ=np.r_[0,np.cumsum((v[1:]+v[:-1])/2)];ax[1].plot(t,integ,color=col,ls=ls,label='Ederer '+method)
ax[1].set(xlabel='Years since diagnosis',ylabel='Integrated survival gap (months)',ylim=(0,9));ax[1].legend(fontsize=8)
for a,l in zip(ax,['a','b']):a.text(-.15,1.04,l,transform=a.transAxes,fontweight='bold',fontsize=12)
fig.tight_layout(w_pad=2);save(fig,1)
s=pd.read_csv(A/'standardized_period_trends.csv');fig,axs=plt.subplots(2,2,figsize=(7.2,5.4));glist=['2.1 Meningiomas','2.2 Tumors of the cranial and paraspinal nerves','Positive histology','Radiography without microscopic confirm'];titles=['Meningiomas','Cranial/paraspinal nerve tumors','Positive histology','Radiography without microscopic confirmation']
for ax,g,title,l in zip(axs.flat,glist,titles,'abcd'):
 for method,col,off in [('I',blue,-.03),('II',orange,.03)]:
  z=s[(s.method==method)&(s.group==g)&s.supported].sort_values('period');xx=np.arange(3)+off;ax.errorbar(xx,z.gap,yerr=[z.gap-z.lo,z.hi-z.gap],fmt='o-',color=col,ms=3,capsize=2,lw=1,label='Ederer '+method)
 ax.set_xticks(range(3),['2004–08','2009–13','2014–18'],fontsize=8);ax.set_title(title,fontsize=8.5);ax.set_ylabel('Five-year gap (months)');ax.text(-.18,1.08,l,transform=ax.transAxes,fontweight='bold',fontsize=12)
axs[0,0].legend(fontsize=7);fig.tight_layout(h_pad=2,w_pad=2);save(fig,3)

plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','axes.linewidth':.7,'legend.frameon':False})
blue='#24658B';gray='#525B64'
def save(fig,n):
 fig.savefig(O/f'Figure_{n}.tif',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'});fig.savefig(O/f'Figure_{n}.svg',bbox_inches='tight');fig.savefig(W/f'Focused_Figure_{n}.png',dpi=160,bbox_inches='tight');plt.close(fig)
groups=['2.1 Meningiomas','2.2 Tumors of the cranial and paraspinal nerves','2.3 Ependymal tumors','2.5 Neuronal and mixed neuronal-glial tumors','2.6 Mesenchymal, nonmeningothelial tumors','C75.1-Pituitary gland','C75.2-Craniopharyngeal duct'];labels=['Meningiomas','Cranial/paraspinal nerve tumors','Ependymal tumors','Neuronal/mixed neuronal-glial tumors','Mesenchymal nonmeningothelial tumors','Pituitary gland C75.1','Craniopharyngeal duct C75.2']
z=pd.read_csv(A/'core_standardized_rmst.csv');z=z[z.method=='I'].set_index('group').loc[groups];y=np.arange(len(groups))
fig,ax=plt.subplots(1,2,figsize=(7.4,4.2),gridspec_kw={'width_ratios':[1.05,1]})
for i,r in enumerate(z.itertuples()):ax[0].plot([r.observed_rmst,r.expected_rmst],[i,i],color='#B4BEC6',lw=1.2)
ax[0].scatter(z.observed_rmst,y,color=gray,s=20,label='Observed');ax[0].scatter(z.expected_rmst,y,color=blue,s=20,marker='s',label='Expected')
ax[0].set(yticks=y,yticklabels=labels,xlabel='Ten-year RMST (months)',xlim=(85,116),ylim=(6.6,-.6));ax[0].legend(loc='lower left',bbox_to_anchor=(0,-.40),ncol=2,fontsize=8)
ax[1].errorbar(z.gap,y,xerr=[z.gap-z.lo,z.hi-z.gap],fmt='o',color=blue,ms=4,capsize=2);ax[1].axvline(0,color=gray,lw=.7,ls=':');ax[1].set(yticks=y,yticklabels=[],xlabel='Expected − observed RMST (months)',xlim=(-3,27),ylim=(6.6,-.6))
for a,l in zip(ax,['a','b']):a.text(-.12,1.04,l,transform=a.transAxes,fontweight='bold',fontsize=11)
fig.tight_layout(w_pad=2);save(fig,1)
p=pd.read_csv(A/'pooled_age_profiles.csv');periods=['0-19','20-39','40-59','60-79','80-84','85-89','90+'];glist=[groups[0],groups[1],groups[5],groups[6]];titles=['Meningiomas','Cranial/paraspinal nerve tumors','Pituitary gland C75.1','Craniopharyngeal duct C75.2'];fig,axs=plt.subplots(2,2,figsize=(7.4,5.7))
for ax,g,title,l in zip(axs.flat,glist,titles,'abcd'):
 z=p[(p.group==g)&p.age_group.isin(periods)].set_index('age_group').reindex(periods);x=np.arange(7);valid=z.supported.fillna(False).to_numpy(bool)
 ax.axvspan(.5,3.5,color='#EAEFF2',zorder=0);ax.axhline(0,color=gray,lw=.7,ls=':');ax.plot(x[valid],z.loc[valid,'gap'],color=blue,lw=.8)
 for i in np.flatnonzero(valid):
  r=z.iloc[i];small=r.n<20;ax.errorbar(i,r.gap,yerr=[[r.gap-r.lo],[r.hi-r.gap]],fmt='o',mfc='white' if small else blue,mec=gray if small else blue,ecolor=gray if small else blue,ms=4,capsize=2,lw=.9)
 for i in np.flatnonzero(~valid):ax.text(i,.03,'NA',transform=ax.get_xaxis_transform(),ha='center',fontsize=7,color=gray)
 ax.set_xticks(x,['0–19','20–39','40–59','60–79','80–84','85–89','90+'],rotation=35,ha='right',fontsize=7);ax.set_title(title,fontsize=9);ax.set_ylabel('Ten-year RMST difference (months)',fontsize=8);ax.text(-.18,1.08,l,transform=ax.transAxes,fontweight='bold',fontsize=11)
fig.tight_layout(h_pad=2,w_pad=2);save(fig,2)
