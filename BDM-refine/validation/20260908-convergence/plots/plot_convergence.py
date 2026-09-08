"""Lightweight vector figures from verified, separately computed measurements."""
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from house_style import PALETTE, use_house_style
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common import ROOT, now, sha, write_json

LABELS={
    'bound_mass':r'$\Delta M_{\rm bound}\ [\%]$',
    'aperture_total_mass':r'$\Delta M_{\rm aperture}\ [\%]$',
    'aperture_radius':r'$\Delta R_{\rm aperture}\ [\%]$',
    'vmax':r'$\Delta V_{\max}\ [\%]$',
    'axis_ba':r'$\Delta(b/a)\ [\%]$',
    'axis_ca':r'$\Delta(c/a)\ [\%]$',
    'bulk_velocity_km_s':r'$|\Delta\boldsymbol{v}|\ [\mathrm{km\,s^{-1}}]$',
    'centre_distance_mpc_h':r'$|\Delta\boldsymbol{x}|\ [h^{-1}\mathrm{Mpc}]$'}
KINDS={'particle':'Particle resolution','force':'Force resolution','time':'Timestep resolution'}


def values(sequence):return np.array(sequence,dtype=float)


def points(ax,x,y,color,label=None,error=None,marker='o'):
    ax.plot(x,y,color=color,lw=1.0,marker=marker,ms=3.2,mec='k',mew=.35,label=label)
    if error is not None:ax.errorbar(x,y,yerr=error,fmt='none',ecolor=color,elinewidth=.7,capsize=0)


def statistic(ax,x,row,key,color,minimum):
    bins=row['matched_statistics'][key]
    counts=values([r['count'] for r in bins]);q=values([r['q16_median_q84'] for r in bins])
    valid=np.array(row['valid_mass_bins'])&(counts>=minimum)
    q[~valid]=np.nan
    ax.fill_between(x,q[:,0],q[:,2],color=color,alpha=.13,lw=0)
    points(ax,x,q[:,1],color)


def figure(rows,edges,kind,z,floor,minimum,supplement=False,partial=False):
    fig,axes=plt.subplots(2,3,figsize=(11.8,6.3),sharex=True)
    axes=axes.ravel();x=(edges[:-1]+edges[1:])/2
    keys=(['aperture_total_mass','axis_ba','bulk_velocity_km_s','centre_distance_mpc_h','unresolved','left_fraction']
          if supplement else ['abundance','bound_mass','vmax','aperture_radius','axis_ca','completeness'])
    for ax,key in zip(axes,keys):
        if key in LABELS:
            ax.set_ylabel(LABELS[key])
            if key not in ['bulk_velocity_km_s','centre_distance_mpc_h']:
                band=2 if key in ['vmax','aperture_radius'] else 5
                ax.axhspan(-band,band,color='.9',zorder=0)
                ax.axhline(0,color='k',lw=.7)
            else:ax.set_ylim(bottom=0)
        elif key=='abundance':
            ax.set_ylabel(r'$100(n_{\rm left}/n_{\rm ref}-1)\ [\%]$')
            ax.axhspan(-5,5,color='.9',zorder=0);ax.axhline(0,color='k',lw=.7)
        elif key=='unresolved':ax.set_ylabel(r'Unresolved $V_{\max}$ [\%]');ax.set_ylim(-2,102)
        else:
            ax.set_ylabel('Reference completeness' if key=='completeness' else 'Left matched fraction')
            ax.axhline(.9,color='.45',lw=.8,ls='--');ax.set_ylim(-.02,1.03)
        ax.grid(True,ls=':',alpha=.25)
        ax.set_xlim(max(12.35,edges[0]),edges[-1])
    handles=[];labels=[]
    for number,row in enumerate(rows):
        color=PALETTE[number%len(PALETTE)];valid=np.array(row['valid_mass_bins'])
        for ax,key in zip(axes,keys):
            if key in LABELS:statistic(ax,x,row,key,color,minimum)
            elif key=='abundance':
                good=valid&(values(row['left_counts'])>=minimum)&(values(row['right_counts'])>=minimum)
                y=100*(values(row['abundance_ratio'])-1);error=100*values(row['abundance_ratio_jackknife8_sigma'])
                y[~good]=np.nan;error[~good]=np.nan
                points(ax,x,y,color,error=error)
            elif key=='unresolved':
                count=values([b['matched_count'] for b in row['vmax_resolution']])
                unresolved=values([b['either_unresolved_count'] for b in row['vmax_resolution']])
                y=np.divide(100*unresolved,count,out=np.full_like(count,np.nan),where=count>0)
                y[~(valid&(count>=minimum))]=np.nan;points(ax,x,y,color)
            else:
                field='reference_completeness' if key=='completeness' else 'left_matched_fraction'
                count='reference_eligible_counts' if key=='completeness' else 'left_eligible_counts'
                y=values(row[field]);y[~(valid&(values(row[count])>=minimum))]=np.nan
                points(ax,x,y,color)
        handle,=axes[0].plot([],[],color=color,marker='o',ms=3.2,mec='k',mew=.35)
        handles.append(handle);labels.append(f"{row['coarse']}/{row['reference']}")
    for ax in axes[3:]:ax.set_xlabel(r'$\log_{10}(M_{\rm bound}/[h^{-1}M_\odot])$')
    fig.legend(handles,labels,loc='upper right',bbox_to_anchor=(.97,1.0),ncol=3,fontsize=11)
    context=rf'{KINDS[kind]}, $z={z}$; $N_{{\rm bound}}\geq {floor}$ in both matched haloes'
    if partial:context='INCOMPLETE CAMPAIGN: '+context
    fig.text(.085,.965,context,fontsize=11,ha='left')
    fig.text(.085,.025,
        f'Only entire mass bins above both resolution cuts; at least {minimum} objects per plotted statistic. '
        'Shifts are left/reference minus one.\n'
        'Shading: matched 16--84 percentile scatter; abundance errors: paired eight-octant jackknife. '
        'Grey bands are reference scales, not certified accuracy.',fontsize=8.2)
    fig.subplots_adjust(left=.085,right=.98,bottom=.14,top=.9,wspace=.33,hspace=.12)
    return fig


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=ROOT/'convergence.json')
    parser.add_argument('--floor',type=int,choices=[100,300,1000],default=300)
    parser.add_argument('--minimum-count',type=int,default=30)
    args=parser.parse_args();data=json.loads(args.input.read_text())
    use_house_style();plt.rcParams.update({'font.size':10,'axes.labelsize':11,'xtick.labelsize':9,'ytick.labelsize':9})
    output=ROOT/'figures'/f'bdm_convergence_n{args.floor}.pdf';output.parent.mkdir(exist_ok=True)
    selected=[r for r in data['comparisons'] if r['particle_floor']==args.floor]
    if not selected:raise RuntimeError('No measured comparison is available to plot')
    pages=[]
    staged=output.with_name(output.stem+'.tmp.pdf')
    with PdfPages(staged) as pdf:
        for z in [2,1,0]:
            for kind in KINDS:
                rows=[r for r in selected if r['redshift']==z and r['kind']==kind]
                if not rows:continue
                for supplement in [False,True]:
                    fig=figure(rows,values(data['log10_mass_edges']),kind,z,args.floor,args.minimum_count,
                               supplement=supplement,partial=not data['completed'])
                    pdf.savefig(fig);plt.close(fig)
                    pages.append(dict(page=len(pages)+1,redshift=z,kind=kind,supplement=supplement,
                                      pairs=[r['coarse']+'/'+r['reference'] for r in rows]))
    staged.replace(output)
    report=dict(completed=True,complete_campaign=data['completed'],created_at_utc=now(),
                input=str(args.input),input_sha256=sha(args.input),script_sha256=sha(__file__),
                particle_floor=args.floor,minimum_count=args.minimum_count,
                pdf=str(output),pdf_sha256=sha(output),pages=pages,
                interpretation='Reference bands illustrate 5% mass/shape/abundance and 2% radius/Vmax scales; they do not establish physical accuracy.')
    write_json(ROOT/f'plot-manifest-n{args.floor}.json',report)
    print(output,len(pages),'pages')


if __name__=='__main__':main()
