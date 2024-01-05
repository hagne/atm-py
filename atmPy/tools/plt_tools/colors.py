import colorsys as _colorsys
import numpy as _np
import matplotlib.pylab as _plt


def wavelength_to_rgb(wavelength, gamma=0.8, rgbrange=1):
    '''This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''

    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0

    R *= rgbrange
    G *= rgbrange
    B *= rgbrange
    #     return (int(R), int(G), int(B))

    return (R, G, B)



def _hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    rgb = _np.array([int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3)]) / 255
    return rgb


def _rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


class Color(object):
    def __init__(self, color, model='hsv', color_scale = 1):
        """
        Parameters
        ----------
        color: depending ond model
        model: (['hsv'], 'rgb', 'hex')
        color_scale:
            ignored when model == 'hex'
        """
        if type(color).__name__ in ['list', 'tuple']:
            color = _np.array(color)

        if model != 'hex':
            color = color.astype(_np.float)
            color /= color_scale

        if len(color) == 4:
            color, alpha = color[:-1], color[-1]
        else:
            alpha = 1

        if model == 'rgb':
            color = _np.array(_colorsys.rgb_to_hsv(*color))

        elif model == 'hls':
            rgb = _colorsys.hls_to_rgb(*color)
            color = _np.array(_colorsys.rgb_to_hsv(*rgb))

        elif model == 'hex':
            rgb = _hex_to_rgb(color)
            if len(rgb) == 4:
                rgb = rgb[:-1]
            color = _np.array(_colorsys.rgb_to_hsv(*rgb))

        self._hsv = color
        self.alpha = alpha

    @property
    def rgb(self):
        #         self.model = 'rgb'
        rgb = _np.array(_colorsys.hsv_to_rgb(*self.hsv))
        return rgb

    @property
    def rgba(self):
        #         self.model = 'rgb'
        rgb = _np.array(_colorsys.hsv_to_rgb(*self.hsv))
        rgba = _np.zeros(4)
        rgba[:-1] = rgb
        rgba[-1] = self.alpha
        return rgba

    @property
    def brightness(self):
        return self._hsv[2]

    @brightness.setter
    def brightness(self, value):
        assert (0 <= value <= 1)
        self._hsv[2] = value

    @property
    def hsv(self):
        return self._hsv

    @property
    def hue(self):
        return self._hsv[0]

    @hue.setter
    def hue(self, value):
        #         assert(0<= value <= 1)
        self._hsv[0] = value

    @property
    def saturation(self):
        return self._hsv[1]

    @saturation.setter
    def saturation(self, value):
        assert (0 <= value <= 1)
        self._hsv[1] = value

    def get_sitance2other_color(self, other, show_colors = False):
        """
        Returns distance to other color. Not clear what this is based on, but
        works for finding a matching color, that might be slidely differnt.

        Parameters
        ----------
        other : Color instance
            the color to compare to.

        Returns
        -------
        float: color distance

        """
        
        rgb1 = self.rgb
        rgb2 = other.rgb
        rm = 0.5*(rgb1[0]+rgb2[0])
        d = sum((2+rm,4,3-rm)*(rgb1-rgb2)**2)**0.5
        if show_colors:
            self.show(scale = 0.5)
            other.show(scale = 0.5)
        return d

    def show(self, alpha=True, scale = 1):
        f, a = _plt.subplots()
        f.set_size_inches(f.get_size_inches() * scale)
        a.axis('off')
        if alpha:
            f.set_facecolor(self.rgba)
        else:
            f.set_facecolor(self.rgb)
        return a