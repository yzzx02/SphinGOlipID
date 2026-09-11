"""Plot only audited theoretical-library totals, never sample identifications."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullLocator

OUT=Path(__file__).resolve().parent


def main():
    audit=json.loads((OUT/'library_totals.json').read_text(encoding='utf-8'))
    own=audit['SphinGOlipID']
    msdial=audit['MS-DIAL']['unique_candidates']
    lipidin=audit['LipidIN']['unique_candidates']
    lda=audit['LDA2']['unique_species_candidates']
    plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.linewidth':1.1,
                         'svg.fonttype':'none','savefig.facecolor':'white'})
    fig,(a,b)=plt.subplots(1,2,figsize=(14,6.7),gridspec_kw={'width_ratios':[1,1.8]})
    fig.subplots_adjust(left=.075,right=.97,bottom=.23,top=.84,wspace=.6)
    fig.suptitle('Sphingolipid theoretical library size',fontsize=18,fontweight='bold',y=.96)
    fig.text(.5,.905,'Source-specific candidate spaces; experimental identifications are not counted',
             ha='center',fontsize=11,color='#555555')
    vals=[own['manuscript_ms1'],own['file_records'],own['unique_names']]
    bars=a.bar(range(3),vals,color=['#A8B9CA','#7998B7','#507493'],width=.6,
               edgecolor='#48617A',linewidth=.7)
    for bar,v in zip(bars,vals):
        a.text(bar.get_x()+bar.get_width()/2,v+1700,f'{v:,}',ha='center',fontsize=11,fontweight='bold')
    a.set_xticks(range(3),['Manuscript\nreported','Provided CSV\nrecords','Provided CSV\nunique names'])
    a.set_ylim(0,101000)
    a.set_ylabel('Number of MS1 entries')
    a.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{int(v):,}'))
    a.set_title('A   SphinGOlipID MS1 library',loc='left',fontweight='bold',pad=18)
    labels=['SphinGOlipID*','MS-DIAL\nVS69','LipidIN\nlocal archive','LDA2','LipiDetective']
    values=[own['expansion'],msdial,lipidin]
    colors=['#93ACC5','#83B6EF','#75C9C3']
    y=[4,3,2,1,0]
    b.set_xscale('log')
    b.set_xlim(10000,50000000)
    for yy,v,c in zip(y,values,colors):
        b.barh(yy,v-10000,left=10000,height=.58,color=c,edgecolor='#586B7C',linewidth=.7,
               hatch='///' if yy==4 else None)
        b.text(v*1.11,yy,f'{v:,}',va='center',ha='left',fontweight='bold',fontsize=11)
    b.text(15000,1,f'Configurable; example MS1 list = {lda:,}',va='center',fontsize=10,color='#686081')
    b.text(15000,0,'Not applicable: sequence prediction model',va='center',fontsize=10,color='#876451')
    b.set_yticks(y,labels)
    b.set_ylim(-.55,4.65)
    b.set_xlabel('Number of unique theoretical candidates (log scale)',labelpad=10)
    b.xaxis.set_major_locator(LogLocator(base=10,numticks=5))
    b.xaxis.set_minor_locator(NullLocator())
    b.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{int(v):,}'))
    b.set_title('B   Chain combinations and reference candidates',loc='left',fontweight='bold',pad=18)
    for ax in [a,b]:
        ax.spines[['top','right']].set_visible(False)
        ax.tick_params(axis='both',labelsize=10)
    b.spines['left'].set_visible(False)
    b.tick_params(axis='y',length=0,pad=12)
    fig.text(.075,.137,'* SphinGOlipID: raw chain combinations from the supplied CSV, before chemical-feasibility filtering.',fontsize=10)
    fig.text(.075,.103,'LCB C12–C30 / DB 0–5; FA C2–C42 / DB 0–3. Two distinct FA positions use the same range in three-chain families.',fontsize=10)
    fig.text(.075,.069,'Reference libraries: polarity, adduct and hierarchical-level duplicates removed by library candidate identity. Native ranges differ.',fontsize=10)
    fig.text(.075,.035,'LDA2 example lists are not its capacity limit. LipiDetective has no fixed theoretical search-library size. Missing bars do not mean zero.',fontsize=10)
    fig.savefig(OUT/'sphingolipid_theoretical_library_totals.png',dpi=300)
    fig.savefig(OUT/'sphingolipid_theoretical_library_totals.svg')
    plt.close(fig)
    rows=[
        ['SphinGOlipID','MS1 manuscript-reported entries',own['manuscript_ms1'],'Reported; differs from supplied CSV'],
        ['SphinGOlipID','MS1 supplied CSV records',own['file_records'],'Includes duplicate names'],
        ['SphinGOlipID','MS1 supplied CSV unique names',own['unique_names'],'Exact name deduplication only'],
        ['SphinGOlipID','Raw chain combinations',own['expansion'],'User-approved paper ranges; no chemical-feasibility filtering'],
        ['MS-DIAL','VS69 unique sphingolipid reference candidates',msdial,'Positive and negative combined; native library ranges'],
        ['LipidIN','Local archive unique sphingolipid reference candidates',lipidin,'Positive and negative combined; LEVEL and adduct duplicates removed'],
        ['LDA2','Recommended example MS1 list candidates',lda,'Example configuration only; software-wide total is not fixed'],
        ['LipiDetective','Fixed theoretical search library','','Not applicable; sequence prediction model'],
    ]
    with (OUT/'theoretical_library_totals.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(['Software','Metric','Total','Definition']);writer.writerows(rows)


if __name__=='__main__': main()
