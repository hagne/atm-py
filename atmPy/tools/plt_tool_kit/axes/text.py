from string import ascii_lowercase as _ascii_lowercase

def set_subplots_labels(a, pos = (0.06, 0.9), cycle_start = 0, bbox = {'fc':[1,1,1,0.7], 'ec':'black', 'boxstyle':'round'}):
    """Adds a letter to each plot in a.

    Parameters
    ----------
    pos: tuple
        position of text in each plot
    cycle_start: int
        If 0 the first plot will get (a), if 1 (b), and so on.
    bbox: dict
        parameters for the box around the text.

    Returns
    -------
    list of text instances
    list of box instances
    """
    if type(a).__name__ == 'AxesSubplot':
        a = [a]
    abc = ['({})'.format(c) for c in _ascii_lowercase]#['(a)', '(b)','(c)', '(b)','(c)', '(b)','(c)', '(b)','(c)', '(b)','(c)', '(b)','(c)']
    texts = []
    boxes = []
    for e, at in enumerate(a):
        if at == None:
            continue
        txtat = at.text(pos[0], pos[1], abc[e + cycle_start], fontsize='large', transform = at.transAxes, bbox = bbox)
        bbp=txtat.get_bbox_patch()
        texts.append(txtat)
        boxes.append(bbp)
    return texts, boxes