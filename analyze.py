# Statistical analysis of nonmalignant CNS tumor survival
# Software version 2.0.0
import shutil
from pathlib import Path
import os,csv,re,json
import numpy as np,pandas as pd
W=Path(os.environ.get('PRIVATE_WORK_DIR','work/private')).resolve();W.mkdir(parents=True,exist_ok=True)
O=Path(os.environ.get('ANALYSIS_OUTPUT_DIR','outputs')).resolve();A=O/'Analysis';A.mkdir(parents=True,exist_ok=True)
SRC=Path(os.environ['SEER_EXPORT_DIR']); EXTRA=Path(os.environ.get('SEER_SURVIVAL_EXPORT_DIR',str(SRC/'survival_exports')))

AGE='Age recode with <1 year olds and 90+'; REC='SEER Brain and CNS Recode'; YEAR='Year of diagnosis'; CONF='Diagnostic Confirmation'
def formats(p):
 out={};k=None
 for l in p.read_text(encoding='utf-8-sig').splitlines():
  if l.startswith('[Format='):k=l[8:-1];out[k]={}
  elif l.startswith('['):k=None
  elif k and '=' in l:
   a,b=l.split('=',1);out[k][a]=b.strip('"').strip()
 return out

def read(p):
 fm=formats(p.with_suffix('.dic'));head=None;rows=[]
 for r in csv.reader(p.open(encoding='utf-8-sig'),delimiter='\t'):
  if r and r[0]=='Page type':head=[s.replace('~',' ') for s in r]
  elif head and 'Interval' in head and len(r)==len(head) and r[0] in ('0','Life Page','Survival Life Page'):rows.append(dict(zip(head,r)))
 d=pd.DataFrame(rows);cols=list(d.columns[1:d.columns.get_loc('Interval')])
 for c in cols:d[c]=d[c].map(lambda x:fm.get(c,{}).get(x,x).strip())
 def interval(s):
  if s.isdigit():return int(s)
  if s in ('Time 0','0 mo'):return 0
  if re.fullmatch(r'<\s*1 mo',s):return 1
  m=re.search(r'<\s*(\d+)',s)
  return int(m[1]) if m else np.nan
 d['Interval']=d.Interval.map(interval)
 for c in d.columns:
  if c not in cols+['Page type']:d[c]=pd.to_numeric(d[c],errors='coerce')
 return d,cols

def calc(d,tau):
 d=d[d.Interval<=tau].sort_values('Interval');t=d.Interval.to_numpy();s=d['Observed Survival (Cum)'].to_numpy();e=d['Expected Survival (Cum)'].to_numpy()
 if not np.array_equal(t,np.arange(tau+1)) or not np.isfinite(s+e).all():return None
 n=d['Adjusted Alive'].to_numpy();dead=d.Died.to_numpy();g=np.cumsum(np.divide(dead,n*(n-dead),out=np.zeros_like(dead,dtype=float),where=(n>dead)&(n>0)))
 w=np.ones(tau+1);w[[0,-1]]=.5;cov=np.outer(s,s)*g[np.minimum.outer(np.arange(tau+1),np.arange(tau+1))];se=np.sqrt(w@cov@w);gap=w@(e-s)
 return dict(n=int(d.iloc[0]['Alive at Start']),tau=tau,gap=gap,se=se,lo=gap-1.96*se,hi=gap+1.96*se,os=s[-1],expected=e[-1],rmst=w@s,expected_rmst=w@e)

cache={}
def load(n,original=False):
 key=(n,original)
 if key not in cache:cache[key]=read((SRC if original else EXTRA)/(n+'.txt'))
 return cache[key]
records=[];cells=[];curves=[]
files=[('I','overall10','Nonmalignant_CNS_RS_EdererI_10y_dx2004_2013_overall',120),('I','overall5','Nonmalignant_CNS_RS_EdererI_5y_dx2004_2018_overall_Life',60),('I','overall15','Nonmalignant_CNS_RS_EdererI_15y_dx2004_2008_overall',180),('II','overall10','nonmalignant_CNS_RS_10y_overall',120),('II','overall5','nonmalignant_CNS_RS_5y_overall',60),('II','overall15','nonmalignant_CNS_RS_15y_overall',180)]
for method,domain,name,tau in files:
 d,cols=load(name,method=='II');r=calc(d,tau);records.append(dict(method=method,domain=domain,**r));curves.append(d.assign(method=method,domain=domain))
 if domain=='overall10':
  for h in [12,60]:records.append(dict(method=method,domain='primary'+str(h),**calc(d,h)))
for method in ['I','II']:
 for domain,name in [('recode','nonmalignant_CNS_RS_10y_brain_recode_age5_Ederer'+method),('site','Nonmalignant_CNS_RS_Ederer'+method+'_10y_C75_site_age5'),('year_recode','Nonmalignant_CNS_RS_Ederer'+method+'_5y_dxyear_brain_recode_age5'),('year_confirmation','Nonmalignant_CNS_RS_Ederer'+method+'_5y_dxyear_confirmation_age5')]:
  d,cols=load(name);tau=60 if 'year' in domain else 120
  for k,z in d.groupby(cols,sort=False):
   k=k if isinstance(k,tuple) else (k,);labels=dict(zip(cols,k));n=z.loc[z.Interval==0,'Alive at Start'].sum()
   if n<=0:continue
   r=calc(z,tau);cells.append(dict(method=method,domain=domain,**labels,n=n,complete=r is not None,**({k:v for k,v in r.items() if k!='n'} if r else {})))
  print('parsed',method,domain,len(d),flush=True)
ct=pd.DataFrame(cells);ct.to_csv(A/'age_cell_estimates.csv',index=False);pd.DataFrame(records).to_csv(A/'overall_estimates.csv',index=False)
# Explicit common support 20-79 years; no interpolation through unsupported tails.
main=pd.read_csv(EXTRA/'Nonmalignant_CNS_10y_cohort_audit.txt',sep='\t',dtype=str)
# Exact primary-cohort members and normalized record keys.
for col in ['Patient ID','Record number recode']:
 main[col]=main[col].str.replace(',','',regex=False).str.lstrip('0')
assert len(main)==main['Patient ID'].nunique()==212236
main['age']=pd.to_numeric(main['Age recode with single ages and 90+'].str.extract(r'(\d+)',expand=False));age0=main.age
weight=[]
for lo in range(20,80,5):weight.append((f'{lo:02d}-{lo+4:02d} years',int(age0.between(lo,lo+4).sum())))
weights=pd.DataFrame(weight,columns=['age','n']);weights['weight']=weights.n/weights.n.sum();weights.to_csv(A/'standard_age_weights.csv',index=False);wd=dict(zip(weights.age,weights.weight))
std=[]
for method in ['I','II']:
 for domain,gcol in [('recode',REC),('site','Primary Site - labeled')]:
  q=ct[(ct.method==method)&(ct.domain==domain)]
  if gcol not in q or q[gcol].isna().all():gcol=next(c for c in q.columns if 'Primary Site' in c) if domain=='site' else REC
  for g,z in q.groupby(gcol):
   z=z[z[AGE].isin(wd)].copy();valid=z.complete.all() and set(z[AGE])==set(wd)
   row=dict(method=method,domain=domain,group=g,age_min=20,age_max=79,n=int(z.n.sum()),supported=valid)
   if valid:
    w=z[AGE].map(wd).to_numpy();gap=w@z.gap;se=np.sqrt((w*w)@(z.se*z.se));row.update(gap=gap,se=se,lo=gap-1.96*se,hi=gap+1.96*se)
   std.append(row)
pd.DataFrame(std).to_csv(A/'standardized_tumor_gaps.csv',index=False)
# Period estimates: baseline-count weighted mixture of yearly curves within age, then fixed age weights.
trends=[]
for method in ['I','II']:
 for domain,gcol in [('year_recode',REC),('year_confirmation',CONF)]:
  q=ct[(ct.method==method)&(ct.domain==domain)&ct[YEAR].str.fullmatch(r'20\d\d',na=False)&ct[AGE].isin(wd)].copy();q['year']=q[YEAR].astype(int);q['period']=pd.cut(q.year,[2003,2008,2013,2018],labels=['2004-2008','2009-2013','2014-2018'])
  for (group,period),z in q.groupby([gcol,'period'],observed=True):
   if domain=='year_recode' and group not in ['2.1 Meningiomas','2.2 Tumors of the cranial and paraspinal nerves']:continue
   if domain=='year_confirmation' and group=='Microscopically confirmed':continue # parent summary overlaps children
   age=[]
   for ag,zz in z.groupby(AGE):
    # Require all exported nonempty yearly cells to be complete; no silent dropping.
    if not zz.complete.all():continue
    ww=zz.n/zz.n.sum();age.append(dict(age=ag,gap=ww@zz.gap,var=(ww*ww)@(zz.se*zz.se),n=int(zz.n.sum())))
   a=pd.DataFrame(age);valid=not a.empty and set(a.age)==set(wd)
   row=dict(method=method,domain=domain,group=group,period=str(period),supported=valid,n=int(z.n.sum()))
   if valid:
    ww=a.age.map(wd);gap=ww@a.gap;se=np.sqrt((ww*ww)@a['var']);row.update(gap=gap,se=se,lo=gap-1.96*se,hi=gap+1.96*se)
   trends.append(row)
pd.DataFrame(trends).to_csv(A/'standardized_period_trends.csv',index=False)
# Conditional observed survival: exact session curves, full covariance delta method.
cond=[]
for name in ['nonmalignant_CNS_RS_10y_overall','nonmalignant_CNS_RS_10y_brain_recode','nonmalignant_CNS_RS_10y_C75_by_primary_site']:
 d,cols=load(name,True)
 for key,z in (d.groupby(cols,sort=False) if cols else [('Overall',d)]):
  z=z.sort_values('Interval');r=calc(z,120)
  if r is None:continue
  s=z['Observed Survival (Cum)'].to_numpy()[:121];n=z['Adjusted Alive'].to_numpy()[:121];dead=z.Died.to_numpy()[:121];g=np.cumsum(np.divide(dead,n*(n-dead),out=np.zeros_like(dead,dtype=float),where=(n>dead)&(n>0)));cov=np.outer(s,s)*g[np.minimum.outer(np.arange(121),np.arange(121))]
  J=np.zeros((61,121))
  for j,i in enumerate(range(60,121)):J[j,i]+=1/s[60];J[j,60]-=s[i]/s[60]**2
  cc=J@cov@J.T;w=np.ones(61);w[[0,-1]]=.5;cs=s[60:]/s[60];se=np.sqrt(cc[-1,-1]);rm=w@cs;rs=np.sqrt(w@cc@w)
  cond.append(dict(group=str(key),n=r['n'],os=cs[-1],lo=cs[-1]-1.96*se,hi=cs[-1]+1.96*se,rmst=rm,rmst_lo=rm-1.96*rs,rmst_hi=rm+1.96*rs))
pd.DataFrame(cond).to_csv(A/'conditional_survival.csv',index=False)
# Private intermediate records remain outside the release.
main.to_pickle(W/'private_main.pkl')
summary={'primary_n':len(main),'standard_population_n':int(weights.n.sum())}
(W/'analysis_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')

# Period contrasts, source curves and alternative age ranges.
from pathlib import Path
import csv,json,re
import pandas as pd,numpy as np

s=pd.read_csv(A/'standardized_period_trends.csv');out=[]
for (m,dom,g),z in s[s.supported].groupby(['method','domain','group']):
 if not {'2004-2008','2014-2018'}<=set(z.period):continue
 z=z.set_index('period');d=z.loc['2014-2018','gap']-z.loc['2004-2008','gap'];se=np.hypot(z.loc['2014-2018','se'],z.loc['2004-2008','se']);out.append(dict(method=m,domain=dom,group=g,difference=d,lo=d-1.96*se,hi=d+1.96*se))
pd.DataFrame(out).to_csv(A/'period_contrasts.csv',index=False)
# Export source curves, paired-method integrity and selected same-age sensitivity ranges.
ns={'EXTRA':EXTRA,'SRC':SRC,'read':read,'REC':REC,'AGE':AGE}
rows=[];curve=[]
for method in ['I','II']:
 p=(ns['EXTRA']/('Nonmalignant_CNS_RS_EdererI_10y_dx2004_2013_overall.txt')) if method=='I' else ns['SRC']/'nonmalignant_CNS_RS_10y_overall.txt'
 d,cols=ns['read'](p);curve.append(d[['Interval','Observed Survival (Cum)','Expected Survival (Cum)','Alive at Start','Died','Lost to Followup']].assign(method=method))
pd.concat(curve).to_csv(A/'primary_monthly_curves.csv',index=False)
ct=pd.read_csv(A/'age_cell_estimates.csv');wt=pd.read_csv(A/'standard_age_weights.csv')
for lo,hi in [(20,74),(40,79)]:
 ww=wt[wt.age.str[:2].astype(int).between(lo,hi)].copy();ww.weight=ww.n/ww.n.sum();wd=dict(zip(ww.age,ww.weight))
 for (method,domain,group),z in ct[ct.domain.isin(['recode','site'])].assign(group=ct[ns['REC']].fillna(ct.get('Primary Site - labeled'))).groupby(['method','domain','group']):
  z=z[z[ns['AGE']].isin(wd)]
  if not z.complete.all() or set(z[ns['AGE']])!=set(wd):continue
  w=z[ns['AGE']].map(wd);gap=w@z.gap;se=np.sqrt((w*w)@(z.se*z.se));rows.append(dict(method=method,domain=domain,group=group,age_min=lo,age_max=hi,gap=gap,lo=gap-1.96*se,hi=gap+1.96*se))
pd.DataFrame(rows).to_csv(A/'alternative_age_ranges.csv',index=False)
print(pd.DataFrame(out).query("method=='I'").to_string(index=False))


# Exact-member fixed-reference resampling.
import sys,os

# Read the primary actuarial life table for observed-curve validation.
import csv,pandas as pd,numpy as np,json
from pathlib import Path
m=pd.read_pickle(W/'private_main.pkl');t=pd.to_numeric(m['Number of Intervals (Calculated)']).to_numpy();status=m['End Calc Vital Status (Adjusted)'].to_numpy()
dead=np.bincount(t[status=='Dead'],minlength=121)[1:121];lost=np.bincount(t[status=='Untraced'],minlength=121)[1:121];surv=int((status=='Alive').sum());start=len(m)-np.r_[0,np.cumsum(dead+lost)[:-1]];adj=start-lost/2;s=np.r_[1,np.cumprod(1-dead/adj)]
h=None;rows=[]
for r in csv.reader(open(SRC/'nonmalignant_CNS_RS_10y_overall.txt',encoding='utf-8-sig'),delimiter='\t'):
 if r and r[0]=='Page type':h=r
 elif h and 'Interval' in h and len(r)==len(h) and r[0]=='0':rows.append(dict(zip(h,r)))
x=pd.DataFrame(rows);obs=x['Observed Survival (Cum)'].astype(float).to_numpy();print('max curve difference',np.max(abs(s-obs)))
print('death count maxdiff',np.max(abs(dead-x.Died.astype(int).to_numpy()[1:])),'lost maxdiff',np.max(abs(lost-x['Lost to Followup'].astype(int).to_numpy()[1:])))
if np.max(abs(s-obs))<.000051:
 counts=np.r_[dead,lost,surv];N=counts.sum();rng=np.random.default_rng(20261024);samples=rng.multinomial(N,counts/N,size=2000);dd=samples[:,:120];ll=samples[:,120:240];aa=N-np.c_[np.zeros(2000),np.cumsum(dd+ll,axis=1)[:,:-1]]-ll/2;ss=np.c_[np.ones(2000),np.cumprod(1-dd/aa,axis=1)];w=np.ones(121);w[[0,-1]]=.5
 # fixed-reference resampling, never labelled full pipeline
 E=[];h=None
 for r in csv.reader(open(EXTRA/'Nonmalignant_CNS_RS_EdererI_10y_dx2004_2013_overall.txt',encoding='utf-8-sig'),delimiter='\t'):
  if r and r[0]=='Page type':h=r
  elif h and 'Interval' in h and len(r)==len(h) and r[0]=='0':E.append(float(r[h.index('Expected Survival (Cum)')]))
 vals=(np.array(E)-ss)@w;cs=ss[:,120]/ss[:,60];crm=(ss[:,60:]/ss[:,60,None])@np.r_[.5,np.ones(59),.5]
 out={'replicates':2000,'seed':20261024,'max_observed_reconstruction_error':float(np.max(abs(s-obs))),'gap_mean':vals.mean(),'gap_sd':vals.std(ddof=1),'gap_ci':np.quantile(vals,[.025,.975]).tolist(),'conditional_os_ci':np.quantile(cs,[.025,.975]).tolist(),'conditional_rmst_ci':np.quantile(crm,[.025,.975]).tolist(),'scope':'Resampling observed interval/status categories with expected curve fixed; not full-pipeline resampling'}
 (A/'bootstrap_validation.json').write_text(json.dumps(out,indent=2));print(out)


# Link morphology to the exact retained primary records.
m=pd.read_pickle(W/'private_main.pkl')
fields=['Patient ID','Record number recode','Primary Site','Histologic Type ICD-O-3','ICD-O-3 Hist/behav','Behavior code ICD-O-3','SEER Brain and CNS Recode','Diagnostic Confirmation']
x=pd.read_csv(SRC/'nonmalignant_CNS_case_listing_2004_2023_SES_treatment.txt',sep='\t',usecols=fields,dtype=str)
for z in [m,x]:
 z['Patient ID']=z['Patient ID'].str.replace(',','',regex=False).str.lstrip('0');z['Record number recode']=z['Record number recode'].str.lstrip('0')
y=m.merge(x,on=['Patient ID','Record number recode'],how='left',validate='one_to_one',suffixes=('','_raw'))
assert len(y)==212236 and y['Histologic Type ICD-O-3'].notna().all()
for c in ['Primary Site','Behavior code ICD-O-3','SEER Brain and CNS Recode','Diagnostic Confirmation']:
 assert y[c].eq(y[c+'_raw']).all(),c
y.to_pickle(W/'private_classified.pkl')

ct=pd.read_csv(A/'age_cell_estimates.csv');wt=pd.read_csv(A/'standard_age_weights.csv');wd=dict(zip(wt.age,wt.weight));out=[]
for (method,domain),q in ct[ct.domain.isin(['recode','site'])].groupby(['method','domain']):
 col=REC if domain=='recode' else 'Primary Site - labeled'
 for group,z in q.groupby(col):
  z=z[z[AGE].isin(wd)]
  if set(z[AGE])!=set(wd) or not z.complete.all():continue
  w=z[AGE].map(wd);obs=w@z.rmst;exp=w@z.expected_rmst;gap=exp-obs;se=np.sqrt((w*w)@(z.se*z.se))
  out.append(dict(method=method,domain=domain,group=group,n=int(z.n.sum()),age_range='20-79',observed_rmst=obs,expected_rmst=exp,gap=gap,se=se,lo=gap-1.96*se,hi=gap+1.96*se,proportional_gap_percent=100*gap/exp))
pd.DataFrame(out).to_csv(A/'core_standardized_rmst.csv',index=False)
age=ct[ct.domain.isin(['recode','site'])].copy();age['group']=age[REC].fillna(age['Primary Site - labeled']);age=age[['method','domain','group',AGE,'n','complete','rmst','expected_rmst','gap','lo','hi','se']];age=age.rename(columns={AGE:'age_stratum','rmst':'observed_rmst'});age.to_csv(A/'age_specific_rmst.csv',index=False)
m=pd.read_pickle(W/'private_classified.pkl');rec=[];cross=[];behavior=[]
for site in ['752','753']:
 for scope,q in [('all ages',m[m['Primary Site'].eq(site)]),('20-79',m[m['Primary Site'].eq(site)&m.age.between(20,79)])]:
  for (hist,label,beh,r),z in q.groupby(['Histologic Type ICD-O-3','ICD-O-3 Hist/behav','Behavior code ICD-O-3',REC]):
   rec.append(dict(site='C'+site[:2]+'.'+site[2:],scope=scope,histology=hist,histology_behavior_label=label,behavior=beh,recode=r,n=len(z),percent=100*len(z)/len(q)))
  for (beh,r),z in q.groupby(['Behavior code ICD-O-3',REC]):cross.append(dict(site='C'+site[:2]+'.'+site[2:],scope=scope,behavior=beh,recode=r,n=len(z),percent=100*len(z)/len(q)))
  hist=q['Histologic Type ICD-O-3'];ben=q['Behavior code ICD-O-3'].eq('Benign');matched=hist.isin(['9350','9351','9352']) if site=='752' else hist.eq('9361')
  behavior.append(dict(site='C'+site[:2]+'.'+site[2:],scope=scope,n=len(q),benign=int(ben.sum()),uncertain=int((~ben).sum()),target_histologies='9350-9352' if site=='752' else '9361',target_n=int(matched.sum()),target_percent=100*matched.mean()))
pd.DataFrame(rec).to_csv(A/'site_histology_behavior.csv',index=False);pd.DataFrame(cross).to_csv(A/'site_recode_crosswalk.csv',index=False);pd.DataFrame(behavior).to_csv(A/'site_composition_summary.csv',index=False)
# Exact retained members support observed RMST sensitivity; no synthetic expected curves.
sens=[]
for site,target in [('752',['9350','9351','9352']),('753',['9361'])]:
 q=m[m['Primary Site'].eq(site)];filters={'All site members':np.ones(len(q),dtype=bool),'Histology restricted':q['Histologic Type ICD-O-3'].isin(target).to_numpy(),'Uncertain behavior only':q['Behavior code ICD-O-3'].eq('Borderline malignancy').to_numpy()}
 for label,mask in filters.items():
  for scope,z in [('all ages',q[mask]),('20-79',q[mask&q.age.between(20,79).to_numpy()])]:
   t=pd.to_numeric(z['Number of Intervals (Calculated)']).to_numpy();st=z['End Calc Vital Status (Adjusted)'];dd=np.bincount(t[st.eq('Dead')],minlength=121)[1:121].astype(float);ll=np.bincount(t[st.eq('Untraced')],minlength=121)[1:121];risk=len(z)-np.r_[0,np.cumsum(dd+ll)[:-1]];adj=risk-ll/2;s=np.r_[1,np.cumprod(1-dd/adj)];g=np.r_[0,np.cumsum(np.divide(dd,adj*(adj-dd),out=np.zeros_like(dd),where=(adj>dd)&(adj>0)))];cov=np.outer(s,s)*g[np.minimum.outer(np.arange(121),np.arange(121))];w=np.ones(121);w[[0,-1]]=.5;rm=w@s;se=np.sqrt(w@cov@w)
   sens.append(dict(site='C'+site[:2]+'.'+site[2:],scope=scope,restriction=label,n=len(z),observed_rmst=rm,lo=rm-1.96*se,hi=rm+1.96*se,expected_rmst='not estimated for restricted subset',gap='not estimated for restricted subset'))
pd.DataFrame(sens).to_csv(A/'site_observed_rmst_sensitivity.csv',index=False)
# Risk counts at five completed years (start of interval 61).
conditions=pd.read_csv(A/'conditional_survival.csv');atrisk=[]
for name in ['nonmalignant_CNS_RS_10y_overall','nonmalignant_CNS_RS_10y_brain_recode','nonmalignant_CNS_RS_10y_C75_by_primary_site']:
 d,cols=read(SRC/(name+'.txt'))
 for key,z in (d.groupby(cols,sort=False) if cols else [('Overall',d)]):
  v=z[z.Interval==61]
  if len(v)==1:atrisk.append(dict(group=str(key),at_risk_after_60_months=int(v.iloc[0]['Alive at Start'])))
conditions.merge(pd.DataFrame(atrisk),on='group',how='left',validate='one_to_one').to_csv(A/'conditional_survival_with_risk.csv',index=False)
# Direct all-age I and II exports, with no reconstruction from incomplete age strata.
crude=[]
patterns=[('recode','I',EXTRA,'*brain_recode*allage*EdererI.txt'),('site','I',EXTRA,'*C75*allage*EdererI.txt'),('recode','II',SRC,'nonmalignant_CNS_RS_10y_brain_recode.txt'),('site','II',SRC,'nonmalignant_CNS_RS_10y_C75_by_primary_site.txt')]
found=[]
for dom,met,root,pattern in patterns:
 paths=list(root.glob(pattern))
 if not paths:continue
 p=paths[-1];found.append(p.name);d,cols=read(p)
 for key,z in d.groupby(cols,sort=False):
  if len(cols)!=1:raise ValueError('All-age export must have one grouping variable')
  label=key[0] if isinstance(key,tuple) else key;r=calc(z,120)
  if r and r['n']>0:crude.append(dict(method=met,domain=dom,group=label,**r))
pd.DataFrame(crude).to_csv(A/'all_age_direct_rmst.csv',index=False)
report={'direct_all_age_files_found':found,'composition':behavior,'primary_reference':'Internal primary-cohort age distribution, not an external standard population','restricted_sensitivity':'Observed RMST only; no new expected trajectories assumed'}
(W/'focus_summary.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2));print(pd.DataFrame(out).query("method=='I'").to_string(index=False));print(pd.DataFrame(sens).to_string(index=False))

domains=[('recode','nonmalignant_CNS_RS_10y_brain_recode_age5_EdererI.txt',REC),('site','Nonmalignant_CNS_RS_EdererI_10y_C75_site_age5.txt','Primary Site - labeled')]
profiles=[];verification=[]
for domain,name,gc in domains:
 d,cols=read(EXTRA/name)
 for group,g in d.groupby(gc):
  children=[]
  for ag,z in g.groupby(AGE):
   n=z.loc[z.Interval==0,'Alive at Start'].sum()
   if n>0:children.append((ag,z,n))
  for label,lo,hi in [('0-19',0,19),('20-39',20,39),('40-59',40,59),('60-79',60,79),('80-84',80,84),('85-89',85,89),('90+',90,999),('all ages',0,999)]:
   parts=[(ag,z,n) for ag,z,n in children if lo<=int(ag[:2])<=hi];n=sum(z[2] for z in parts)
   if not n:continue
   rr=[calc(z,120) for ag,z,nn in parts];row=dict(domain=domain,group=group,age_group=label,n=int(n),supported=all(r is not None for r in rr))
   if row['supported']:
    numeric=['Alive at Start','Died','Lost to Followup','Adjusted Alive'];pooled=sum(z.set_index('Interval')[numeric].loc[range(121)] for ag,z,nn in parts);dead=pooled.Died.to_numpy();risk=pooled['Adjusted Alive'].to_numpy();s=np.cumprod(1-np.divide(dead,risk,out=np.zeros_like(dead,dtype=float),where=risk>0));e=sum(nn/n*z.set_index('Interval')['Expected Survival (Cum)'].loc[range(121)].to_numpy() for ag,z,nn in parts);gg=np.cumsum(np.divide(dead,risk*(risk-dead),out=np.zeros_like(dead,dtype=float),where=(risk>dead)&(risk>0)));cov=np.outer(s,s)*gg[np.minimum.outer(np.arange(121),np.arange(121))];w=np.ones(121);w[[0,-1]]=.5;gap=w@(e-s);se=np.sqrt(w@cov@w);row.update(observed_rmst=w@s,expected_rmst=w@e,gap=gap,lo=gap-1.96*se,hi=gap+1.96*se,se=se)
   profiles.append(row)
pd.DataFrame(profiles).to_csv(A/'pooled_age_profiles.csv',index=False)
q=pd.DataFrame(profiles)
print(q[(q.group.isin(['2.1 Meningiomas','2.2 Tumors of the cranial and paraspinal nerves','C75.1-Pituitary gland','C75.2-Craniopharyngeal duct']))&q.age_group.ne('all ages')].to_string(index=False))

from pathlib import Path
import numpy as np
import pandas as pd
import json

def run_stability_audit(O,W,SRC,read,calc,AGE,REC):
 A=O/'Analysis';extra=EXTRA
 wt=pd.read_csv(A/'standard_age_weights.csv');wd=dict(zip(wt.age,wt.weight))
 neuronal='2.5 Neuronal and mixed neuronal-glial tumors'
 groups=['2.1 Meningiomas','2.2 Tumors of the cranial and paraspinal nerves','2.3 Ependymal tumors',neuronal,'2.6 Mesenchymal, nonmeningothelial tumors','C75.1-Pituitary gland','C75.2-Craniopharyngeal duct']
 curves={};con=[];land=[]
 for domain,name,gc in [('recode','nonmalignant_CNS_RS_10y_brain_recode_age5_EdererI.txt',REC),('site','Nonmalignant_CNS_RS_EdererI_10y_C75_site_age5.txt','Primary Site - labeled')]:
  d,_=read(extra/name)
  for (group,ag),z in d.groupby([gc,AGE]):
   if group not in groups or ag not in wd:continue
   z=z.sort_values('Interval').set_index('Interval');r=calc(z.reset_index(),120)
   assert r is not None
   curves[group,ag]=(z,r)
   if group==neuronal:
    con.append(dict(age_stratum=ag,n=r['n'],weight=wd[ag],deaths_0_120=int(z.Died.sum()),at_risk_60=int(z.loc[61,'Alive at Start']),at_risk_96=int(z.loc[97,'Alive at Start']),at_risk_start_last_interval=int(z.loc[120,'Alive at Start']),observed_rmst=r['rmst'],expected_rmst=r['expected_rmst'],gap=r['gap'],lo=r['lo'],hi=r['hi'],weighted_gap_contribution=wd[ag]*r['gap'],weighted_variance=wd[ag]**2*r['se']**2))
 c=pd.DataFrame(con).sort_values('age_stratum');c['gap_contribution_percent']=100*c.weighted_gap_contribution/c.weighted_gap_contribution.sum();c['variance_percent']=100*c.weighted_variance/c.weighted_variance.sum();c.to_csv(A/'neuronal_age_contributions.csv',index=False)
 for group in groups:
  rs=[]
  for ag in wd:
   z,r=curves[group,ag];rs.append((wd[ag],z,r))
  # Ratio of separately standardized observed and expected curves, consistently at both landmarks.
  row=dict(group=group,n=sum(r['n'] for w,z,r in rs),method='I',age_range='20-79')
  for month in [60,120]:
   obs=sum(w*z.loc[month,'Observed Survival (Cum)'] for w,z,r in rs);exp=sum(w*z.loc[month,'Expected Survival (Cum)'] for w,z,r in rs)
   var=0
   for w,z,r in rs:
    zz=z.loc[:month];nn=zz['Adjusted Alive'].to_numpy();dd=zz.Died.to_numpy();gg=np.divide(dd,nn*(nn-dd),out=np.zeros_like(dd,dtype=float),where=(nn>dd)&(nn>0)).sum();var+=w*w*z.loc[month,'Observed Survival (Cum)']**2*gg
   se=np.sqrt(var)/exp;row.update({f'observed_{month}':obs,f'expected_{month}':exp,f'relative_survival_{month}':obs/exp,f'rs_lo_{month}':obs/exp-1.96*se,f'rs_hi_{month}':obs/exp+1.96*se})
  ob=sum(w*r['rmst'] for w,z,r in rs);ex=sum(w*r['expected_rmst'] for w,z,r in rs);se=np.sqrt(sum(w*w*r['se']**2 for w,z,r in rs));row.update(observed_rmst=ob,expected_rmst=ex,gap=ex-ob,lo=ex-ob-1.96*se,hi=ex-ob+1.96*se);land.append(row)
 pd.DataFrame(land).to_csv(A/'relative_survival_rmst_comparison.csv',index=False)
 # Alternative age ranges omitting sparse older cells.
 alt=pd.read_csv(A/'alternative_age_ranges.csv');alt=alt[~alt.age_max.isin([59,69])]
 added=[]
 for lo,hi in [(20,59),(20,69)]:
  ww=wt[wt.age.str[:2].astype(int).between(lo,hi)];weights=dict(zip(ww.age,ww.n/ww.n.sum()))
  for group in groups:
   vals=[(w,curves[group,ag][1]) for ag,w in weights.items()];gap=sum(w*r['gap'] for w,r in vals);se=np.sqrt(sum(w*w*r['se']**2 for w,r in vals));added.append(dict(method='I',domain='site' if group.startswith('C') else 'recode',group=group,age_min=lo,age_max=hi,gap=gap,lo=gap-1.96*se,hi=gap+1.96*se))
 pd.concat([alt,pd.DataFrame(added)],ignore_index=True).to_csv(A/'alternative_age_ranges.csv',index=False)
 m=pd.read_pickle(W/'private_classified.pkl');composition=[]
 for group in groups:
  q=m[m['Primary Site'].eq(group[1:3]+group[4:5])] if group.startswith('C') else m[m[REC].eq(group)]
  scopes=[('all ages',q),('20-79',q[q.age.between(20,79)]),('60-79',q[q.age.between(60,79)])]+[(f'{lo}-{lo+4}',q[q.age.between(lo,lo+4)]) for lo in [60,65,70,75]]
  for scope,z in scopes:
   for (hist,label,behavior),v in z.groupby(['Histologic Type ICD-O-3','ICD-O-3 Hist/behav','Behavior code ICD-O-3']):composition.append(dict(group=group,age_scope=scope,histology=hist,histology_behavior_label=label,behavior=behavior,n=len(v),percent=100*len(v)/len(z)))
 pd.DataFrame(composition).to_csv(A/'major_group_composition.csv',index=False)
 summary={'neuronal_gap':float(c.weighted_gap_contribution.sum()),'age60_79_n':int(c[c.age_stratum.str[:2].astype(int)>=60].n.sum()),'age60_79_gap_contribution':float(c[c.age_stratum.str[:2].astype(int)>=60].weighted_gap_contribution.sum()),'age60_79_variance_percent':float(c[c.age_stratum.str[:2].astype(int)>=60].variance_percent.sum()),'age60_79_weight':float(c[c.age_stratum.str[:2].astype(int)>=60].weight.sum()),'rs_definition':'ratio of separately age-standardized observed and Ederer I expected survival curves; not weighted mean of age-specific ratios'}
 (A/'stability_audit_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8');print(json.dumps(summary,indent=2));print(c.to_string(index=False));print(pd.concat([alt,pd.DataFrame(added)]).query('method=="I" and group==@neuronal').to_string(index=False))



run_stability_audit(O,W,SRC,read,calc,AGE,REC)

from pathlib import Path
import numpy as np
import pandas as pd

def run_nerve_sensitivity(O,W):
 A=O/'Analysis';m=pd.read_pickle(W/'private_classified.pkl');q=m[m['SEER Brain and CNS Recode'].str.startswith('2.2')].copy();wt=pd.read_csv(A/'standard_age_weights.csv');cells=[];summ=[]
 definitions=[('All nerve group members',q),('Positive histology',q[q['Diagnostic Confirmation'].eq('Positive histology')]),('Radiography without microscopic confirmation',q[q['Diagnostic Confirmation'].eq('Radiography without microscopic confirm')]),('Benign behavior /0',q[q['Behavior code ICD-O-3'].eq('Benign')]),('Uncertain behavior /1',q[q['Behavior code ICD-O-3'].eq('Borderline malignancy')])]
 for label,part in definitions:
  for age in wt.itertuples():
   lo=int(age.age[:2]);z=part[part.age.between(lo,lo+4)];n=len(z);row=dict(restriction=label,age_stratum=age.age,n=n,weight=age.weight,supported=False)
   if n:
    t=pd.to_numeric(z['Number of Intervals (Calculated)']).to_numpy();status=z['End Calc Vital Status (Adjusted)'];dd=np.bincount(t[status.eq('Dead')],minlength=121)[1:121].astype(float);ll=np.bincount(t[status.eq('Untraced')],minlength=121)[1:121].astype(float);risk=n-np.r_[0,np.cumsum(dd+ll)[:-1]];adj=risk-ll/2
    row.update(deaths=int(dd.sum()),at_risk_60=int(risk[60]),at_risk_start_last_interval=int(risk[-1]))
    if np.all(adj>0):
     so=np.r_[1,np.cumprod(1-dd/adj)];gg=np.r_[0,np.cumsum(np.divide(dd,adj*(adj-dd),out=np.zeros_like(dd),where=adj>dd))];cov=np.outer(so,so)*gg[np.minimum.outer(np.arange(121),np.arange(121))];v=np.ones(121);v[[0,-1]]=.5;rm=v@so;se=np.sqrt(v@cov@v);row.update(supported=True,observed_rmst=rm,se=se,lo=rm-1.96*se,hi=rm+1.96*se)
   cells.append(row)
  cc=pd.DataFrame([r for r in cells if r['restriction']==label]);valid=bool(cc.supported.all());row=dict(restriction=label,age_range='20-79',n=int(cc.n.sum()),all_age_n=len(part),complete_age_cells=int(cc.supported.sum()),supported=valid,expected_rmst='not estimated for restricted subset',gap='not estimated for restricted subset')
  if valid:
   rm=cc.weight@cc.observed_rmst;se=np.sqrt((cc.weight**2)@(cc.se**2));row.update(observed_rmst=rm,se=se,lo=rm-1.96*se,hi=rm+1.96*se)
  summ.append(row)
 pd.DataFrame(cells).to_csv(A/'nerve_observed_sensitivity_age_cells.csv',index=False);pd.DataFrame(summ).to_csv(A/'nerve_observed_sensitivity.csv',index=False)
 records=pd.concat([pd.DataFrame(cells).assign(record_type='age cell'),pd.DataFrame(summ).assign(record_type='standardized summary')],ignore_index=True);records.to_csv(O/'Supplementary_Data_12_Nerve_observed_sensitivity.csv',index=False)
 # Reconstruction agrees with the reported observed component to the precision of exported curves.
 ref=pd.read_csv(A/'core_standardized_rmst.csv').query('method=="I" and group.str.startswith("2.2")',engine='python').iloc[0]
 assert abs(summ[0]['observed_rmst']-ref.observed_rmst)<0.006
 print(pd.DataFrame(summ).to_string(index=False))


run_nerve_sensitivity(O,W)
from pathlib import Path
import pandas as pd
import numpy as np
import json

def run_nerve_expected_validation(O,W,SRC,read,calc,AGE,REC):
 A=O/'Analysis';extra=EXTRA;name='nonmalignant_CNS_RS_10y_brain_recode_confirmation_age5_EdererI.txt'
 d,cols=read(extra/name);assert set(cols)=={REC,AGE,'Diagnostic Confirmation'},cols
 original,_=read(extra/'nonmalignant_CNS_RS_10y_brain_recode_age5_EdererI.txt')
 group='2.2 Tumors of the cranial and paraspinal nerves';nerve=d[d[REC].eq(group)].copy();reference_group=original[original[REC].eq(group)]
 leaf=nerve[nerve['Diagnostic Confirmation'].ne('Microscopically confirmed')]
 keys=[AGE,'Interval'];countcols=['Alive at Start','Died','Lost to Followup','Adjusted Alive'];pooled=leaf.groupby(keys)[countcols].sum();ref=reference_group.set_index(keys)[countcols]
 pooled=pooled.reindex(ref.index);assert np.allclose(pooled,ref,equal_nan=True)
 # Reconstruct observed curves from additive interval risk/death counts rather than averaging subgroup survival.
 maxdiff=0;expected_diffs=[]
 for ag,z in pooled.reset_index().groupby(AGE):
  z=z.sort_values('Interval');rr=z['Adjusted Alive'].to_numpy();dd=z.Died.to_numpy();so=np.cumprod(1-np.divide(dd,rr,out=np.zeros_like(dd,dtype=float),where=rr>0));oo=reference_group[reference_group[AGE].eq(ag)].sort_values('Interval')['Observed Survival (Cum)'].to_numpy();valid=np.isfinite(oo);maxdiff=max(maxdiff,float(np.max(np.abs(so[valid]-oo[valid]))))
  parts=[]
  for confirm,zz in leaf[leaf[AGE].eq(ag)].groupby('Diagnostic Confirmation'):
   zz=zz.set_index('Interval').reindex(range(121));nn=zz.loc[0,'Alive at Start']
   if nn>0:parts.append((nn,zz['Expected Survival (Cum)'].to_numpy()))
  if parts:
   ep=np.array([v for n,v in parts]);weights=np.array([n for n,v in parts]);full=np.isfinite(ep).all(axis=0);reference=reference_group[reference_group[AGE].eq(ag)].set_index('Interval').reindex(range(121))['Expected Survival (Cum)'].to_numpy();full&=np.isfinite(reference);recombined=weights@ep/weights.sum();expected_diffs.extend(np.abs(recombined[full]-reference[full]).tolist())
 assert maxdiff<0.000051
 assert expected_diffs and max(expected_diffs)<0.000101,'Expected reference does not reproduce original Ederer I baseline mixture'
 wt=pd.read_csv(A/'standard_age_weights.csv');wd=dict(zip(wt.age,wt.weight));m=pd.read_pickle(W/'private_classified.pkl');m=m[m[REC].eq(group)];cells=[];summ=[];monthly=[]
 labels=[('Positive histology','Positive histology'),('Radiography without microscopic confirm','Radiography without microscopic confirmation')]
 for exported,label in labels:
  part=nerve[nerve['Diagnostic Confirmation'].eq(exported)];records=m[m['Diagnostic Confirmation'].eq(exported)];support=True
  assert int(part.loc[part.Interval.eq(0),'Alive at Start'].sum())==len(records)
  for ag,w in wd.items():
   z=part[part[AGE].eq(ag)].sort_values('Interval');n=int(z.loc[z.Interval.eq(0),'Alive at Start'].sum());lo=int(ag[:2]);raw=records[records.age.between(lo,lo+4)];assert len(raw)==n
   t=pd.to_numeric(raw['Number of Intervals (Calculated)']).to_numpy();status=raw['End Calc Vital Status (Adjusted)'];dead=np.bincount(t[status.eq('Dead')],minlength=121)[1:121];lost=np.bincount(t[status.eq('Untraced')],minlength=121)[1:121]
   assert np.array_equal(z[z.Interval.between(1,120)].Died.to_numpy(),dead)
   assert np.array_equal(z[z.Interval.between(1,120)]['Lost to Followup'].to_numpy(),lost)
   r=calc(z,120);row=dict(restriction=label,age_stratum=ag,weight=w,n=n,supported=r is not None)
   if r:row.update(observed_rmst=r['rmst'],expected_rmst=r['expected_rmst'],gap=r['gap'],se=r['se'],lo=r['lo'],hi=r['hi'],deaths=int(dead.sum()),at_risk_60=int(z.loc[z.Interval.eq(61),'Alive at Start'].iloc[0]),at_risk_start_last_interval=int(z.loc[z.Interval.eq(120),'Alive at Start'].iloc[0]))
   else:support=False
   cells.append(row);monthly.append(z.assign(restriction=label))
  cc=pd.DataFrame([r for r in cells if r['restriction']==label]);row=dict(restriction=label,age_range='20-79',n=int(cc.n.sum()),all_age_n=len(records),supported=support,complete_age_cells=int(cc.supported.sum()),method='I')
  if support:
   ob=cc.weight@cc.observed_rmst;ex=cc.weight@cc.expected_rmst;se=np.sqrt((cc.weight**2)@(cc.se**2));row.update(observed_rmst=ob,expected_rmst=ex,gap=ex-ob,se=se,lo=ex-ob-1.96*se,hi=ex-ob+1.96*se,proportional_gap_percent=100*(ex-ob)/ex)
  summ.append(row)
 assert all(r['supported'] for r in summ),'Do not silently omit unsupported age cells'
 core=pd.read_csv(A/'core_standardized_rmst.csv');r=core[(core.method=='I')&core.group.eq(group)].iloc[0];summ.insert(0,dict(restriction='All nerve group members',age_range='20-79',n=int(r.n),all_age_n=len(m),supported=True,complete_age_cells=12,method='I',observed_rmst=r.observed_rmst,expected_rmst=r.expected_rmst,gap=r.gap,se=r.se,lo=r.lo,hi=r.hi,proportional_gap_percent=r.proportional_gap_percent))
 obs=pd.read_csv(A/'nerve_observed_sensitivity.csv')
 for row in summ:
  independent=obs[obs.restriction.eq(row['restriction'])].iloc[0];assert abs(row['observed_rmst']-independent.observed_rmst)<0.006
 pd.DataFrame(cells).to_csv(A/'nerve_confirmation_rmst_age_cells.csv',index=False);pd.DataFrame(summ).to_csv(A/'nerve_confirmation_rmst.csv',index=False)
 pd.concat(monthly,ignore_index=True).to_csv(A/'nerve_confirmation_monthly_curves.csv',index=False)
 pd.concat([pd.DataFrame(cells).assign(record_type='age estimate'),pd.DataFrame(summ).assign(record_type='standardized summary')],ignore_index=True).to_csv(O/'Supplementary_Data_13_Nerve_confirmation_RMST.csv',index=False)
 report={'input_file':name,'nerve_full_age_n':len(m),'pooled_leaf_interval_counts_match':'all counts match original across all exported ages; microscopy parent excluded from sums','observed_curve_max_absolute_reconstruction_difference':maxdiff,'expected_baseline_mixture_comparable_points':len(expected_diffs),'expected_baseline_mixture_max_absolute_difference':max(expected_diffs),'exact_primary_subset_n_deaths_censors_match':'passed in all 24 age cells','same_internal_age_weights':True,'complete_cells':24,'observed_rmst_reconstruction_tolerance_months':0.006,'interpretation':'Subgroup-specific Ederer I references; confirmation restriction changes membership and does not establish causal effects or detection bias. Fixed expected curves and weights for confidence intervals.'}
 (A/'nerve_confirmation_validation.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2));print(pd.DataFrame(summ).to_string(index=False))


run_nerve_expected_validation(O,W,SRC,read,calc,AGE,REC)
pd.read_csv(A/'core_standardized_rmst.csv').to_csv(O/'Supplementary_Data_1_Core_RMST.csv',index=False)
cr=pd.read_csv(A/'all_age_direct_rmst.csv');over=pd.read_csv(A/'overall_estimates.csv').query("domain=='overall10'").assign(group='Overall cohort');pd.concat([cr,over],ignore_index=True).to_csv(O/'Supplementary_Data_2_All_age_RMST.csv',index=False)
fine=pd.read_csv(A/'age_specific_rmst.csv').assign(estimate_scope='fine age stratum');pooled=pd.read_csv(A/'pooled_age_profiles.csv').assign(method='I',estimate_scope='pooled actual age group').rename(columns={'age_group':'age_stratum','supported':'complete'});pd.concat([fine,pooled],ignore_index=True).to_csv(O/'Supplementary_Data_3_Age_specific_RMST.csv',index=False)
shutil.copy2(A/'site_histology_behavior.csv',O/'Supplementary_Data_4_Site_classification.csv');shutil.copy2(A/'standard_age_weights.csv',O/'Supplementary_Data_5_Internal_age_weights.csv');shutil.copy2(A/'conditional_survival_with_risk.csv',O/'Supplementary_Data_6_Conditional_survival.csv')
period=pd.read_csv(A/'standardized_period_trends.csv').assign(record_type='period estimate');delta=pd.read_csv(A/'period_contrasts.csv').rename(columns={'difference':'gap'}).assign(period='2014-2018 minus 2004-2008',record_type='period difference');pd.concat([period,delta],ignore_index=True).to_csv(O/'Supplementary_Data_7_Exploratory_periods.csv',index=False)
shutil.copy2(A/'primary_monthly_curves.csv',O/'Supplementary_Data_8_Overall_figure_source.csv')
shutil.copy2(A/'neuronal_age_contributions.csv',O/'Supplementary_Data_9_Neuronal_stability.csv')
shutil.copy2(A/'relative_survival_rmst_comparison.csv',O/'Supplementary_Data_10_Relative_survival_RMST.csv')
shutil.copy2(A/'major_group_composition.csv',O/'Supplementary_Data_11_Major_group_composition.csv')
