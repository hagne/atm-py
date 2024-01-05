import numpy as _np

def pick_line_and_highlight(a, lw = 2, scale2fit = False):
    """Highlights a line in a plot when selected with mouse click. Arrow up and down allows cycling through lines.

    Notes
    -----
    In jupyter use "%matplotlib nbagg" for this to work.

    Attributes
    ----------
    a: matplotlib.Axes instance
    lw: float
        linewidth enhancement factor"""
    f = a.get_figure()
    for g in a.get_lines():
        g.set_picker(5)

    last_artist = {'g': None}

    def apply_change(artist, last_artist):
        #         global last_artist
        globals()['artist'] = artist
        globals()['last_artist'] = last_artist

        if last_artist['g']:

            la = last_artist['g']
            la.set_linewidth(last_artist['lw'])
            la.set_zorder(last_artist['zo'])

        if artist == last_artist['g']:
            last_artist['g'] = None
        else:
            last_artist['g'] = artist
            last_artist['lw'] = artist.get_linewidth()
            last_artist['zo'] =  artist.get_zorder()

            artist.set_linewidth(artist.get_linewidth() * lw)
            artist.set_zorder(100)
        if a.legend_:
            a.legend()

        if scale2fit:
            x, y = artist.get_xydata().transpose()
            dx = abs(_np.nanmax(x) - _np.nanmin(x))
            dy = abs(_np.nanmax(y) - _np.nanmin(y))
            scale = 0.05
            a.set_xlim(_np.nanmin(x) - (dx * scale), _np.nanmax(x) + (dx * scale))
            a.set_ylim(_np.nanmin(y) - (dy * scale), _np.nanmax(y) + (dy * scale))

    def onclick(event, last_artist):
        globals()['event'] = event
        artist = event.artist
        apply_change(artist, last_artist)


    def on_key(event, last_artist):
        #         global last_artist
        lines = a.get_lines()
        globals()['event'] = event
        #         lines = last_artist.axes.get_lines()
        if last_artist['g']:
            idx = lines.index(last_artist['g'])
        else:
            idx = -1
        globals()['idx_was'] = idx
        if event.key == 'up':
            idx -= 1
        elif event.key == 'down':
            idx += 1
        else:
            return

        if idx < 0:
            idx = len(lines) + idx
        elif idx == len(lines):
            idx = 0
        globals()['idx_is'] = idx
        globals()['line'] = lines[idx]
        apply_change(lines[idx], last_artist)



    cid = f.canvas.mpl_connect('pick_event', lambda x: onclick(x, last_artist))
    cid = f.canvas.mpl_connect('key_press_event', lambda x: on_key(x, last_artist))