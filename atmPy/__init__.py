#from . import aerosols
#from . import data_archives
#from . import tools
#from . import file_io as open_file_by_format
#from .file_io import open_atmpy as open_file

# on anaconda importing basemap can cause errors which can be fixed by setting
# the following variable
try:
    import mpl_toolkits.basemap as __
except KeyError as er:
    if er.args[0] == 'PROJ_LIB':
        import os as _os
        _os.environ['PROJ_LIB']  = _os.environ['CONDA_PREFIX']
        import mpl_toolkits.basemap as __
    else:
        raise