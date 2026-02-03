# optional_imports.py
class OptionalImport:
    def __init__(self, name, submodules = None):
        self.module_available = False
        self.module = None
        self.name = name

        self.submodules = submodules
        
        self._attempt_import()
        self._attempt_import_submods()

    def _attempt_import_submods(self):
        if (not isinstance(self.submodules, type(None))) and self.module_available:
            submodules = self.submodules
            
            if not isinstance(submodules, list):
                submodules = [submodules,]
                
            for mod in submodules:
                __import__(f'{self.name}.{mod}')
            

    def _attempt_import(self):
        try:
            self.module = __import__(self.name)
            self.module_available = True
        except ImportError:
            self.module_available = False

    def __getattr__(self, item):
        if not self.module_available:
            raise ImportError(f"{self.name} is required for this feature. Please install it to use this functionality.")
        return getattr(self.module, item)

# Creating the pandas facade
statsmodels = OptionalImport('statsmodels', submodules = ['api','nonparametric.smoothers_lowess'])

#Todo: remove those and replace with the submodule kwarg
# statsmodels_api = OptionalImport('statsmodels.api')
# statsmodels_nonparametric_smoothers_lowess = OptionalImport('statsmodels.nonparametric.smoothers_lowess')
# statsmodels_robust =  OptionalImport('statsmodels.robust')

timezonefinder = OptionalImport('timezonefinder')


mpl_toolkits_basemap = OptionalImport('mpl_toolkits.basemap')
metpy = OptionalImport('metpy')
pptx = OptionalImport('pptx')

ephem = OptionalImport('ephem')

pysolar = OptionalImport('pysolar')

geopy = OptionalImport('geopy', submodules = 'distance')

cartopy = OptionalImport('cartopy', submodules = 'io.img_tiles')
pygam = OptionalImport('pygam')

numba = OptionalImport('numba')
pvlib = OptionalImport('pvlib')

matplotlib = OptionalImport('matplotlib', submodules='pyplot')

pytz = OptionalImport('pytz')

psutil = OptionalImport('psuitl')