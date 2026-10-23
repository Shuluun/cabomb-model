"""Figures for the model-generated waiting-time comparison."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def plot_bridge(output_dir: str | Path) -> Path:
    """Plot each design on its own scale, with units and scenario labels."""
    output=Path(output_dir)
    data=pd.read_csv(output/'bridge_curves.csv')
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True)
    names={'unmodified':'Unmodified','ocn':'OCN','ysa':'YSA','ocn-ysa':'OCN + YSA'}
    for ax,(design,name) in zip(axes.flat,names.items()):
        for beta in (.5,1.,2.):
            frame=data[(data.design==design)&(data.beta==beta)]
            ax.plot(frame.time_h,frame.gross_response,label=f'Shape {beta:g}')
        ax.set_title(name)
        ax.set_xlabel('Time since injection (h)')
        ax.set_ylabel('Response (normalized carrier equivalents)')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
        ax.legend(fontsize=8)
    fig.suptitle('Release-waiting shape at a shared 4 h mean\n200 nm; medium EPR/MPS; fixed V1.3 internalization input')
    fig.tight_layout()
    path=output/'bridge_response.png'
    fig.savefig(path,dpi=170)
    plt.close(fig)
    (output/'figure_caption.txt').write_text(
        'Each panel uses its own y-axis scale. Curves are assumed cumulative waiting responses '
        'convolved with model-generated internalization; units are normalized carrier equivalents, '
        'not measured cargo or free calcium. Source dynamics and mean waiting time are held fixed.\n',encoding='utf-8')
    return path
