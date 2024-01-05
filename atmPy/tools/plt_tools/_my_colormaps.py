import plt_tools
import numpy as _np

def clouds():
    bmax = plt_tools.colors.Color((0 ,0 ,0.1))
    bmin = plt_tools.colors.Color((0 ,0 ,1))
    colors = _np.array([bmin.rgb , bmax.rgb])
    cmap = plt_tools.colormap.creat_cmap(colors=colors, norm = 'linear',
                                         log_min= 0.1,
                                         reverse=False)
    cmap.set_bad((1 ,1 ,1))
    return cmap

def temperature(limit_to = None):
    """use slice() for limit_to"""
    brutal_heat = plt_tools.colors.Color((191 ,19 ,0), color_scale=255, model='rgb')
    #     hot = plt_tools.colors.Color((255,112,3), color_scale= 255, model='rgb')
    warm = plt_tools.colors.Color((255 ,255 ,92), color_scale= 255, model='rgb')
    moderate = plt_tools.colors.Color((41 ,153 ,41), color_scale= 255, model='rgb')
    cold = plt_tools.colors.Color((23 ,16 ,158), color_scale= 255, model='rgb')
    freezing = plt_tools.colors.Color((246 ,207 ,255), color_scale= 255, model='rgb')

    colors = _np.array([freezing.rgb, cold.rgb, moderate.rgb, warm.rgb, brutal_heat.rgb])
    if limit_to:
        colors = colors[limit_to]
    cmap = plt_tools.colormap.creat_cmap(colors=colors, norm = 'linear',
                                         log_min= 0.1,
                                         reverse=False)
    cmap.set_bad((1 ,1 ,1))
    return cmap

def relative_humidity(limit_to = None, reverse = False):
    """use slice() for limit_to"""
    bone_dry = plt_tools.colors.Color((255 ,255 ,196), color_scale=255, model='rgb')
    dry = plt_tools.colors.Color((227 ,200 ,82), color_scale= 255, model='rgb')
    nice = plt_tools.colors.Color((191 ,93 ,20), color_scale= 255, model='rgb')
    wet = plt_tools.colors.Color((41 ,153 ,41), color_scale= 255, model='rgb')
    wet.brightness = 0.5
    wet.saturation = 0.7
    dripping = plt_tools.colors.Color((23 ,16 ,158), color_scale= 255, model='rgb')
    #     freezing = plt_tools.colors.Color((246,207,255), color_scale= 255, model='rgb')

    colors = _np.array([bone_dry.rgb,
                        #                        dry.rgb,
                       nice.rgb, wet.rgb, dripping.rgb])
    if limit_to:
        colors = colors[limit_to]
    cmap = plt_tools.colormap.creat_cmap(colors=colors, norm = 'linear',
                                         log_min= 0.1,
                                         reverse=reverse)
    cmap.set_bad((1 ,1 ,1))
    return cmap

def particle_concentration(limit_to = None):
    """use slice() for limit_to"""
    venus = plt_tools.colors.Color((138 ,59 ,0), color_scale=255, model='rgb')
    venus.brightness *= 0.1
    unhealthy = plt_tools.colors.Color((191 ,19 ,0), color_scale=255, model='rgb')
    #     hot = plt_tools.colors.Color((255,112,3), color_scale= 255, model='rgb')
    moderate = plt_tools.colors.Color((255 ,255 ,92), color_scale= 255, model='rgb')
    #     moderate = plt_tools.colors.Color((41,153,41), color_scale= 255, model='rgb')
    #     cold = plt_tools.colors.Color((23,16,158), color_scale= 255, model='rgb')
    clean = plt_tools.colors.Color((246 ,207 ,255), color_scale= 255, model='rgb')
    clean.hue -= 0.2

    colors = _np.array([clean.rgb, moderate.rgb, unhealthy.rgb, venus.rgb])
    if limit_to:
        colors = colors[limit_to]
    cmap = plt_tools.colormap.creat_cmap(colors=colors, norm = 'linear',
                                         log_min= 0.1,
                                         reverse=False)
    cmap.set_bad((1 ,1 ,1))
    return cmap