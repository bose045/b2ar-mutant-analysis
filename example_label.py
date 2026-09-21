"""Mark every demonstration figure without changing the plotted data."""
_example=False

def set_example(value):
    global _example
    _example=bool(value)

def annotate_example(fig):
    if _example and not getattr(fig,'_example_banner',False):
        fig.text(.5,.995,'SYNTHETIC EXAMPLE — NOT SIMULATION RESULTS',ha='center',va='top',
                 fontsize=10,color='#a02020',bbox={'facecolor':'white','edgecolor':'none','alpha':.95},zorder=1000)
        fig._example_banner=True
