# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/htelg/.spyder2/.temp.py
"""

#ToDo
#- ich denke nicht, dass wir neutral nehen muessen .. unserer laser is polarisiert ...
#- indes of refraction at 405 nm for psl
#- results plotten

# Check using http://omlc.ogi.edu/calc/mie_calc.html


import os
import sys
import types

import matplotlib.cm as mplcm
import matplotlib.colors as colors
import numpy as np
import pylab as plt
from scipy.interpolate import interp1d

from atmPy.aerosols.instruments.POPS import tools
from atmPy.radiation.mie_scattering import bhmie
import pandas as _pd


###########################

def makeMie_diameter(diameterRangeInMikroMeter = [0.1,3.],
            noOfdiameters = 200,
            noOfAngles = 100, # number of scatternig angles
            POPSdesign = 'POPS 2',
            IOR = 1.5,
            WavelengthInUm = .405,
            geometry = "perpendicular",
            #the following was added 20140904
            mirrorJetDist = 10.,
            scale = 'log',
            #below: added 20141030
            broadened = False,
            doreturn = 'tuple'
            ):
    """
    Performs mie_scattering calculations as a function of particle radius

    Arguments
    ---------
    geometry: "perpendicular",

    broadened: {"style": 'gauss',
                 'center': 0.405,
                 'fwhm': 0.005,
                 'cutoff': 2, number of fwhm
                 'noOfwl': 10}
               {'style': 'custom',
                 'spectrum': (wl,intens),
                 'interpolate': 100}
    doreturn: str (['tuple'], 'dict')


    Returns
    -------
    diameter
    array and an array of the intensity of the light scattered onto the detector (this is also 
    not fully correct, since we do not do a decent integration but a sum, at some point it would be nice to change that)
    
    parameters:

    """
    radiusRangeInMikroMeter = np.array(diameterRangeInMikroMeter) / 2
    if scale == 'log':
        dRange = np.logspace(np.log10(radiusRangeInMikroMeter[0]),np.log10(radiusRangeInMikroMeter[1]),noOfdiameters) #radius range 
    elif scale == 'linear':
        dRange = np.linspace(radiusRangeInMikroMeter[0],radiusRangeInMikroMeter[1],noOfdiameters) #radius range 
    elif scale == 'log_20':
        dRange = np.logspace(np.log10(radiusRangeInMikroMeter[0]),np.log10(radiusRangeInMikroMeter[1]),noOfdiameters,base = 100) #radius range 
        
    event = Mie(silent = True, design = POPSdesign, indexOfRef = IOR, diameter = 'dynamic') # Ref. for indexOfRef.  Patterson 2004
    event.set_nang(noOfAngles)
    event.POPSdimensions['mirror(top)-jet distance (mm)'] = float(mirrorJetDist)
    
    singleLine = False
    
    if isinstance(WavelengthInUm,float):
        exWavelengthInUm=np.array([WavelengthInUm])
    elif isinstance(WavelengthInUm,types.ListType):
        exWavelengthInUm=np.array(WavelengthInUm)
    else:
        exWavelengthInUm = WavelengthInUm
    

    if broadened:
        if broadened['style'] == 'gauss':
            wc = broadened['center'] 
            fwhm = broadened['fwhm']
            noOfwl = broadened['noOfwl']
            cutoff = broadened['cutoff']
            # fwxm = wc / ()
            # fwxm = sigmas*np.sqrt(sigmas*np.log(2))* sigma /(sigmas*np.sqrt(sigmas*np.log(10)))
            exWavelengthInUm = np.linspace(wc - (cutoff * fwhm), wc+ (cutoff * fwhm), noOfwl)
            normalizer = tools.gauss_function(exWavelengthInUm, wc, fwhm)
        if broadened['style'] == 'custom':
            exWavelengthInUm = broadened['spectrum'][0]
            normalizer = broadened['spectrum'][1]
            if broadened['interpolate']:
                spectrum = interp1d(exWavelengthInUm,normalizer)
                noOfwl = broadened['interpolate']
                wl_new = np.linspace(exWavelengthInUm.min(),exWavelengthInUm.max(),noOfwl)
                exWavelengthInUm = wl_new
                normalizer = spectrum(wl_new)
        lightsource_spectrum = _pd.DataFrame(normalizer, index = exWavelengthInUm)
                
    
    else:
        normalizer = np.ones(exWavelengthInUm.shape)/float(len(exWavelengthInUm))
        noOfwl = len(exWavelengthInUm)
        
    if len(exWavelengthInUm) == 1:
        singleLine = True
        
    output = np.zeros((exWavelengthInUm.shape[0]+1,dRange.shape[0]))
    for e,i in enumerate(exWavelengthInUm):
        event.set_wavelength(i)
        perpInt = []
        for i in dRange:
            event.set_r(i)
            intensity =  event.get_detectableIntensity(geometry)
            perpInt.append(intensity)
        diameter = np.array(2 * np.array(dRange))
        scatteringIntensity = np.array(perpInt)
#         if broadened:     
        output[0]+= normalizer[e] * scatteringIntensity/normalizer.sum()
        if len(exWavelengthInUm) == 1:
            return diameter, scatteringIntensity
        elif not singleLine: #why am I doing that?
            output[e+1]=scatteringIntensity
        elif e == int(len(exWavelengthInUm)/2): #why am I doing that?
            output[1] = scatteringIntensity
        
    if doreturn == 'dict':
        out = {'diameter': diameter,
               'scattering_int': output[0],
               'lightsource_spectrum': lightsource_spectrum}
        return out #diameter,output#,(exWavelengthInUm,normalizer)
    elif doreturn == 'tuple':
        return diameter, output
    
###########################################################    
class Mie():
    """ Creates a Mie object
    An introduction to mie_scattering scattering: http://omlc.org/education/ece532/class3/mie_math.html
    
    What to do next:
    1. Define the following parameters:
    \t silent:    \t True or False - if True a lot of output will be generated, 
                  \t usefull for debugging or understanding what whas actually done
    \t diameter:  \t number or array type - unit is micro meter - if array, 
                  \t scattering efficiency of all diameters in array will be calculated. 
                  \t Note, diameter and wavelength can not be both array type
    \t wavelength:\t number or array type - unit is micro meter - if array, 
                  \t scattering efficiency of all wavelength in array will be calculated. 
                  \t Note, diameter and wavelength can not be both array type.
    \t indexOfRef:\t number or array type - if array, scattering efficiency for all IORs in array will be calculated. 
                  \t ToDo, sortout diameter wavelength array issues
    \t design:    \t string - this sets all parameters related to the geometry of the particular version of POPS. 
                  \t see docstring of set_dimensions for details
    \t nang:      \t integer - angle steps in which calculations are performed.
    \t material:  \t Not fully functional. TODO: create database
    """
    
    def __init__(self, silent = True, diameter = False, wavelength = False, material = False, indexOfRef = False, design = 'POPS 2', nang = 100):
        # if indexOfRef and material:
#             raise TypeError('index of refraction is calculated on the defined material. Therefore, either indexOfRef or material has to be set to False') 
            
        self.material = material
        self.set_wavelength(wavelength)
        self.set_d(diameter)
        self.set_n(n = indexOfRef)
        self.set_nang(nang)
               
        self.silent = silent
        
        self.POPSdimensions = {}
        self.design = design
        if design:
            self.set_dimensions(design)
            
        self.upToDate = False
        self.upToDate_geometry = False
        self.dontAskAgain = False
            
    def set_dimensions(self, design):
        if design == 'POPS 1': # This is the prototype, printed with the old 3D Printer
            mDiameter = self.POPSdimensions['mirror diameter (mm)'] = 25. 
            rMirror = 20. # radius of the underlying sphere of the spherical mirror
            self.POPSdimensions['mirror(top)-jet distance (mm)'] = rMirror - tools.segment_hight(rMirror, mDiameter)
            self.POPSdimensions['angle: jet-mirrorNormal (rad)'] = np.pi/2. # This is the mean scattering angle (90 deg)

        elif design == 'POPS 2': # This is the prototype, printed with the old 3D Printer 
            mDiameter = self.POPSdimensions['mirror diameter (mm)'] = 25. 
            rMirror = 20. # radius of the underlying sphere of the spherical mirror
            self.POPSdimensions['mirror(top)-jet distance (mm)'] = 7.68
            self.POPSdimensions['angle: jet-mirrorNormal (rad)'] = np.pi/2. # This is the mean scattering angle (90 deg)
            self.POPSdimensions['polarization (laser)'] = 'perpendicular'
            
        if not self.silent:
            print("You are using %s" % design)
            print("self.POPSdimensions['mirror(top)-jet distance (mm)']",
                  self.POPSdimensions['mirror(top)-jet distance (mm)'])
            print("self.POPSdimensions['mirror diameter (mm)']", self.POPSdimensions['mirror diameter (mm)'])
            print("self.POPSdimensions['angle: jet-mirrorNormal (rad)']",
                  self.POPSdimensions['angle: jet-mirrorNormal (rad)'])
            print(
            "(not used yet) self.POPSdimensions['polarization (laser)']", self.POPSdimensions['polarization (laser)'] )

            print('Distance from top to bottom of mirror', tools.segment_hight(rMirror, mDiameter))
        
        if self.POPSdimensions['angle: jet-mirrorNormal (rad)'] != np.pi/2:
            raw_input("angle is not 90 deg, you better check everythin, it might not work right")
        self.upToDate_geometry = False
        
    def set_r(self, r):
        if not r:
            print(r)
            self.r = 0.5 * float(raw_input('please define particle diameter in um: '))
        elif np.all(r > 10) and not self.dontAskAgain:
            antwort = raw_input("""This is probably an error, are you sure you are giving the particle radius in um?!? (y)""")
            if antwort.lower() == 'n':
                sys.exit('as you wish ... exit')
            else:
                self.dontAskAgain = True
        else:
            self.r = r
        self.upToDate = False
        
    def set_d(self, d):
        if type(d) == str:
            self.r = d
        else:
            self.set_r(d/2.)
        
    def set_wavelength(self, wavelength):
        if wavelength > 10.:
            antwort = raw_input("""This is probably an error, are you sure you are giving the wavelength in um?!?""")
            if antwort.lower() == 'n':
                sys.exit('as you wish ... exit!')
        self.wavelength = wavelength
        self.upToDate = False
        
    def set_n(self, n = False):
        if n and self.material:
            raise TypeError('Scattering event is set to "material", n can therefore not be changed')
#        elif n and self.n:
#            raise TypeError('Refractive index was defined at the begining, n can therefore not be changed')
        elif self.material:
            materialList = []
            materialList.append({'name' : 'polystyrene', 'get' :  tools.refIndex_polystyrene})
            found = False
            for i in materialList:
                if i['name'] == self.material:
                    self.n = i['get'](self.wavelength)
                    if found:
                        raise ValueError("Found two matching materials! Not possible!")
                    else:
                        found = True
            if not found:
                raise KeyError('%s is not in the list of available materials' % self.material)
        elif not n:
            self.n = complex(raw_input('please define refractive index: '))
        elif n:
            self.n = n   
        self.upToDate = False
          
    def set_nang(self,nang):
        """set the number of angles between 0 and 90 degries"""
        self.nang = nang
        self.upToDate = False
        
    def set_SizeParameter(self): #r, wavelength):
        self.x = 2*np.pi/self.wavelength * self.r
        return self.x
        
    def check_param_Defined(self):
        if not self.r:
            self.set_r(float(raw_input('please define particle radius in um: ')))
#            antwort = F
        if not self.wavelength:
            self.set_wavelength(float(raw_input('please define wavelength in nm: ')))
#            antwort = False
        if not self.n:
            if self.material:
                pass
            else:
                self.set_n(complex(raw_input('please define defractive index: ')))
#            antwort = False
        if not self.nang:
            self.set_nang(raw_input('please set number of angles between 0 and 90 degrees: '))

#     def update(self):
#         self.check_param_Defined()
#         if not self.upToDate:
#             self.set_SizeParameter()
#             self.do_bhmie()
#             self.upToDate = True
#             
#             #self.set_xAxes()
#             if not self.silent:
#                 print 'updated'
#         else:
#             if not self.silent:
#                 print 'no update needed'
            
    def update_hagen(self):
        self.check_param_Defined()
        if not self.upToDate:
            if self.material:
                self.set_n()
            self.set_SizeParameter()
            self.do_bhmie_hagen()
#             self.calc_Natural()
#             self.calc_Perpendicular()
#             self.calc_Parallel()
            self.set_xAxis()
            self.upToDate = True
            
            #self.set_xAxes()
            if not self.silent:
                print('updated')
        else:
            if not self.silent:
                print('no update needed')
                
    def update_geometry(self):
        if not self.upToDate_geometry:
            self.get_mirror_grid()
            self.upToDate = True
        
                   
    def set_xAxis(self):
        noOfPts = 2 * ((self.nang * 2) - 1)
        self.xAxis = np.linspace(0,2 * np.pi,noOfPts)
        

    
    def print_current_parameter(self):
        """prints all relevant parameters"""
        print('particle diameter: \t\t', 2 * self.r)
        print('particle size parameter: \t', self.x)
        print('wavelength: \t\t\t', self.wavelength)
        print('refractive index: \t\t', self.n)
        print('design:  \t\t\t', self.design)
              
#     def do_bhmie(self):
#         self.s1,self.s2,self.qext,self.qsca,self.qback,self.gsca = bhmie.bhmie(self.x, self.n, self.nang)
        
    def do_bhmie_hagen(self):
        if not self.silent:
            self.print_current_parameter()
        bhh = bhmie.bhmie_hagen(self.x, self.n, self.nang)
        s1,s2,self.qext,self.qsca,self.qback,self.gsca = bhh.return_Values()
#         data = (abs(self.s1))**2#/(np.pi * self.x**2 * self.qsca)
        s1_Reverse = s1[::-1]
        self.s1 = np.concatenate((s1,s1_Reverse)) 
        s2_Reverse = s2[::-1]
        self.s2 = np.concatenate((s2,s2_Reverse))


    def print_Data(self):
        for i in self.xAxis:
            print(i, ' , ', self.YNatural[i])
            
    def get_mirror_grid(self):
        # np.set_printoptions(threshold=np.nan)
        # np.set_printoptions(precision=2)
        
        
        dm = self.POPSdimensions['mirror diameter (mm)']
        h = self.POPSdimensions['mirror(top)-jet distance (mm)']
#         print 'h', h
#         print 'dm', dm
        rSphere = tools.sphereRadius_fromGeometry(h, dm)         # 1.
        alphMax = tools.alphamax_fromGeometry(h, dm)
        sSphere = tools.arc_length(rSphere, h)
#         sSphereA = tools.arc_length_alpha(rSphere, alphMax * 2.)          #for test purposes
#         angleRangeArray, dataRangeArray = tools.find_angleRange(self.POPSdimensions['angle: jet-mirrorNormal (rad)'],alphMax,self.xAxis, data)
        angleRangeArray, angleIndexArray = tools.find_angleRange(self.POPSdimensions['angle: jet-mirrorNormal (rad)'], alphMax, self.xAxis)
        stepWidth = sSphere/len(angleRangeArray)
        
        self.angleIndexArray = angleIndexArray
        
        nn = len(angleIndexArray)
        indexMatrix = np.ones((nn,nn), dtype = int) * angleIndexArray
        angleMatrix = np.ones((nn,nn), dtype = int) * angleRangeArray
        
#         
#         offAngleMatrix = np.empty((nn,nn))
#         offAngleMatrix[:] = 0 #np.NAN
#         
        yArcLenghtMatrix = np.empty((nn,nn))
        yArcLenghtMatrix[:] = 0
#         print 'nn', nn
#         print 'stepWidth', stepWidth
        ArcLengthArray = abs(np.array(list(range(int(-nn / 2), 0, 1)) + list(range(0, int(nn / 2), 1))))
        ArcLengthMatrix = np.ones((nn,nn))
        ArcLengthMatrix = (ArcLengthMatrix * ArcLengthArray).transpose() * stepWidth
        
#         print "ArcLengthMatrix"
#         raw_input(ArcLengthMatrix.astype(int))
        
        
        
        rs = tools.sphereSegment_radius(rSphere, angleMatrix - np.pi / 2.)
#         print h
        ss = tools.arc_length(rs, h)
        
#         print "ss"
#         raw_input(ss.astype(int))
        
        ArcLengthMatrix[ArcLengthMatrix > ss/2.] = np.nan
#         print "ArcLengthMatrix"
#         raw_input(ArcLengthMatrix.astype(int))
        
        self.offAngleMatrix = .5 * tools.segment_angle(rs, ArcLengthMatrix)
        


                	
            
    def get_detectableIntensity(self, polarization = "perpendicular"):
        """ In this function I want to calculate a solid angle which is defined by the mirror and then all the light which is scattered into that angle.
            Parameters:
            \t geometry: polarization of the laser with respect to the plane defined by laser beam and the collection direction. values: perpendiular, paralllel or natural
            The following calculations might seam weard. However, since there is no way of calculating the scattering efficiency analytically as a function of the angle I create a...
            1. claculate the distance from the particle to the edge of the mirror: 
                rSphere = tools.sphereRadius_fromGeometry(h,dm)
                    h: shortest distance from particle to plane defined by the edge of the mirror
                    dm: mirror diameter
            2. calculate the maximum colectable scattering angle alphMax. alphMax is the angle between vector(particle-center of morror) and vector(particle-edge of morror)
                alphMax = tools.alphamax_fromGeometry(h,dm)
            3. calculate the arc length which is cut out of the imaginary circle with radius rSphere by the mirror 
                sSphere = tools.arc_length(rSphere, h)
            4. get the section of the mie_scattering_data based on the angle between incident light and mirror normal (probably 90 deg) and the arc length
                angleRangeArray, dataRangeArray = tools.find_angleRange(self.POPSdimensions['angle: jet-mirrorNormal (rad)'],alphMax,self.xAxis, data)
            5. since we are are looking at the scattering in 3D we also need to go to the side. We will do that in a way that the angle stays the same. 
                What will change is the orientation of the scattering plane and ,therefore ,the polariazation of the light with respect to that plane. 
                Therefor we will go along arcs perpendicular to the arc we calculated before (sSphere) using the same stepwith which is defined by the 
                number of data points along the sSphere:
                stepWidth = sSphere/len(angleRangeArray)
            """
        
        #if isinstance(self, collections.Iterable) and not isinstance(a, types.StringTypes):
        #    print 'bla'
        self.update_hagen()
        self.update_geometry()
#         mirror_grid = self.get_mirror_grid()
#         raw_input('watewatte')
        
        whatList = ('natural', 'parallel', 'perpendicular')
        if polarization not in whatList:
            raise ValueError('Geometry has to be one of the following: "%s", "%s", or "%s"? %s is not an option' % (
            whatList[0], whatList[1], whatList[2], polarization))
        fIdx = self.angleIndexArray[0]
        lIdx = self.angleIndexArray[-1]
        s1Selection = self.s1[fIdx:lIdx+1]
        s2Selection = self.s2[fIdx:lIdx+1]
        
        nn = len(self.angleIndexArray)
        
        s1Matrix = (abs(np.ones((nn,nn)) * s1Selection))**2
        s2Matrix = (abs(np.ones((nn,nn)) * s2Selection))**2
        naturalMatrix = (.5 * s1Matrix) + (.5 * s2Matrix)
#         print 's1Matrix'
#         raw_input(s1Matrix)
        
        if len(s1Selection) != nn:
            raise ValueError('not possible %s %s'%(len(s1Selection),len(self.angleIndexArray)))

        
        if polarization == "parallel":
            IntMatrix = ((np.cos(self.offAngleMatrix))**2 * s2Matrix) + ((np.sin(self.offAngleMatrix))**2 * s1Matrix)  
        
        elif polarization == "perpendicular":
            IntMatrix = ((np.sin(self.offAngleMatrix))**2 * s2Matrix) + ((np.cos(self.offAngleMatrix))**2 * s1Matrix)  
        
        elif polarization == "natural":
            naturalMatrix[np.isnan(self.offAngleMatrix)] = 0
            IntMatrix = naturalMatrix 
        
        IntMatrix[np.isnan(IntMatrix)] = 0
        integratedIntensity = IntMatrix.sum() 

        return integratedIntensity# * stepWidth**2

def plot_polar(dataList, log = False):


    NUM_COLORS = len(dataList)

    cm = plt.get_cmap('gist_rainbow')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
#     ax = fig.add_subplot(111)
    fig, ax = plt.subplots(figsize=(12,9), subplot_kw=dict(projection='polar'))
    # old way:
    #ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    # new way:
    ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
    col = ['r','b']
    for e,i in enumerate(dataList):
        #color = cm(1.*e/NUM_COLORS)  # color will now be an RGBA tuple
        if log:
            ax.plot(i[0], np.log10(i[1]), linewidth=(3 - (2 * e)), color = col[e])
        else:
            ax.plot(i[0], i[1], linewidth=(3 - (2 * e)), color = col[e])
#        ax.plot(i.xAxis, i.YNatural, linewidth=3-e)
#    ax.set_rmax(0.2)
    
    ax.grid(True)		
#    ax.set_title("A line plot on a polar axis", va='bottom')
    if log:
        ax.set_ylim((-3, 0))
    plt.show()

def plot_POPS_calib(dataList, log = (0,0), title=False):
    fig, ax = plt.subplots(figsize=(12,9))
    
    NUM_COLORS = len(dataList)
    cm = plt.get_cmap('gist_rainbow')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
     
    for e,i in enumerate(dataList):
        color = cm(1.*e/NUM_COLORS)  # color will now be an RGBA tuple
        if i['label'] == 'gaus. broadened':
            linewidth = 3
            color = 'black'
        else:
            linewidth = 1
        ax.plot(i['x'], i['y'], label = i['label'], linewidth = linewidth, color = color)
#        ax.plot(i.xAxis, i.YNatural, linewidth=3-e)
#    ax.set_rmax(0.2)
    
    ax.grid(True)	
    if log[1]:	
        ax.set_yscale('log')
    if log[0]:
        ax.set_xscale('log')

#    ax.set_title("A line plot on a polar axis", va='bottom')
#     ax.set_ylim((-3, 0))
    if title:
        ax.set_title(title)
    ax.set_xlabel('Diameter ($\mu$m)')
    ax.set_ylabel('Scattering efficiency (arb. u.)')
    ax.legend()
    plt.show()
    return fig, ax
    
def save(dataList, name ='test', nameAddOn = 'label',extension = '.dat'):
    if not os.path.exists('output'):
        os.makedirs('output')
        print('new directory with name "output" created')
        
    for e,i in enumerate(dataList):
        finalName = 'output/'+name
        finalName += i['label']
        finalName += str(len(i['x'])) + 'pts'
        finalName += extension
        np.savetxt(finalName,np.array([i['x'],i['y']]).transpose())
        if len(finalName) > 35:
            print("Warning!!! filename is longer than 35 characters and therefore might cause trouble with Iogr!")
        
    
##############################################################################################
##############################################################################################
##############################################################################################
 
def default():
    event = Mie()
    event.set_r(.05)
    event.set_wavelength(.405)
    event.set_n(1.5)
    event.set_nang(100)
    event.calc_Natural()

    print('qext: ', event.qext)
    print('qsca: ', event.qsca)
    print('qback: ', event.qback)
    
    plot_polar([(event.xAxis, event.YNatural)], log = False)

def steptest():
    """Result: the programm is actually calculating the radiance -> absolut values are indipendent of number of caluclated points"""
    event = Mie()
    event.set_r(.500)
    event.set_wavelength(.6328)
    event.set_n(1.5 + 0.1j)
    event.set_nang(10)
    event.calc_Natural()

    print('qext: ', event.qext)
    print('qsca: ', event.qsca)
    print('qback: ', event.qback)
    print('natural: ', event.natural)
    
    eventII = Mie()
    eventII.set_r(.500)
    eventII.set_wavelength(.6328)
    eventII.set_n(1.5 + 0.1j)
    eventII.set_nang(50)
    eventII.calc_Natural()

    print('qext: ', event.qext)
    print('qsca: ', event.qsca)
    print('qback: ', event.qback)
    
    plot_polar([event,eventII])

def bhmie_hagen_test():
    event = Mie()
    eventII = Mie()
    r = 1.
    lamb = .405
    n = 1.5 + 0.1j
    noOfPts = 100
    event.set_r(r)
    event.set_wavelength(lamb)
    event.set_n(n)
    event.set_nang(noOfPts)
    event.calc_Natural()

    print('qext: ', event.qext)
    print('qsca: ', event.qsca)
    print('qback: ', event.qback)
    print('natural: ', event.natural)
    
    eventII = Mie()
    eventII.set_r(r)
    eventII.set_wavelength(lamb)
    eventII.set_n(n)
    eventII.set_nang(noOfPts)
    eventII.calc_Natural_hagen()

    print('qext: ', eventII.qext)
    print('qsca: ', eventII.qsca)
    print('qback: ', eventII.qback)
    print('natural: ', eventII.natural)
    
    plot_polar([event,eventII])
    
def test_comparison_to_internet():
    """ this does a comparison of our calculations to that from the online calculator"""
    f = open('data/mieCalculaterOutput.txt','r')
    
    x_ref = []
    y_ref_nat = []
    for line in f.readlines():
        if 'radius' in line:
            print(line)
        elif 'n_real' in line:
            print(line)
        elif 'size parameter' in line:
            print(line)
        elif 'wavelength' in line:
            print(line)
        elif line[0] != '#':
            values = line.split()
            x_ref.append(np.deg2rad(float(values[0])))
            y_ref_nat.append(float(values[1]))
    
    x_ref = np.array(x_ref)  
    y_ref_nat = np.array(y_ref_nat)
    #y_calc = np.array(event.YNatural)
    #sys.exit('ende gut alles gut')    
    event = Mie(wavelength = 0.6328, diameter = 1.0, indexOfRef = 1.5, nang = 100)
    event.calc_Natural_hagen()
    print("length", len(event.xAxis), len(x_ref))
    print(" first x", event.xAxis[25], x_ref[0])
    print('first y', event.YNatural[0], y_ref_nat[0])
    print('max', event.YNatural.max(), y_ref_nat.max())
    plot_polar([(x_ref,y_ref_nat/y_ref_nat.max()),(event.xAxis, event.YNatural/event.YNatural.max())], log = True)
    #plot_polar([(x_ref,y_ref_nat/y_ref_nat.max())])
    
def test_calc_nintyOnly():
    """ plots scattering scattering efficiency at exactly 90 deg compared to the area colected by the mirror of the POPS"""
    
    #event_point = Mie(silent = True, wavelength = 0.6328, diameter = 1.0, indexOfRef = 1.5, nang = 100)
    event_area  = Mie(silent = True, design = 'POPS 1',wavelength = 0.6328, diameter = 1.0, indexOfRef = 1.5, nang = 100)
    
#     event.POPSdimentions['mirror diameter (mm)'] = 25
#     event.POPSdimentions['mirror(top)-jet distance (mm)'] = 12.25
#     event.POPSdimentions['angle: jet-mirrorNormal (rad)'] = np.pi/2.
    
    mDiameter = event_area.POPSdimensions['mirror diameter (mm)'] = 25. 
    rMirror = 20. # radius of the underlying sphere of the spherical mirror
    event_area.POPSdimensions['mirror(top)-jet distance (mm)'] = 1000.#rMirror - tools.segment_hight(rMirror,mDiameter)
    event_area.POPSdimensions['angle: jet-mirrorNormal (rad)'] = np.pi/2. # This is the mean scattering angle (90 deg)
    
    diameters = np.logspace(-1,1,10)
    #xList_point = []
    yList_point = []
    xList_area = []
    yList_area = []
    for i in diameters:
#         print 'i', i
        #event_point.set_d(i)
        #event_point.calc_Natural_hagen()
        #intensity = event_point.YNatural[int(event_point.nang)]
        #xList_point.append(i)
        #yList_point.append(intensity)
        
        event_area.set_d(i)
        event_area.calc_Natural_hagen()
        intensity = event_area.YNatural[int(event_area.nang)]
        yList_point.append(intensity)
        intensity = event_area.get_detectableIntensity()
        xList_area.append(i)
        yList_area.append(intensity)
        
        
    spect_dic_point= {'x': xList_area, 'y': yList_point, 'label': '90'}    
    spect_dic_area = {'x': xList_area, 'y': yList_area, 'label': 'area'}   
    plot_POPS_calib([spect_dic_point, spect_dic_area], log =(1,1))   
    #plot_POPS_calib([spect_dic_point], log =(1,1))
     
def calc_scatteringPattern():
    event = Mie(silent = True)
    # event.set_r(1.00)
    event.set_wavelength(.405)
    event.set_n(1.5)
    event.set_nang(1000)
    
    
    event.POPSdimentions['mirror diameter (mm)'] = 25
    event.POPSdimentions['mirror(top)-jet distance (mm)'] = 12.25
    event.POPSdimentions['angle: jet-mirrorNormal (rad)'] = np.pi/2.
    
    dataList = []
    for i in np.linspace(0.05,0.75,15):        
#         print 'i', i
        event.set_r(i)
        event.calc_Natural_hagen()
        scattPat = event.YNatural.copy()
        xA = event.xAxis.copy()
        dataList.append((xA,scattPat))
        
    plot_polar(dataList, log = True)
    
def calc_intensityAsFktOfRadius_versusRuShan():
    """ creates calibration plot for comparison with Ru-Shans calculations
        -  the aim is to get a pattern at 90 degrees and almost no area of the mirror
    """
    event = Mie(silent = True)
    # event.set_r(1.00)
    event.set_wavelength(.405)
    event.set_n(1.52)
    event.set_nang(1000)
    
    
    event.POPSdimentions['mirror diameter (mm)'] = 25.4
    event.POPSdimentions['mirror(top)-jet distance (mm)'] = 10
    event.POPSdimentions['angle: jet-mirrorNormal (rad)'] = np.pi/2.
    
    xList = []
    yList = []
    for i in np.linspace(0.05,1.5,5000):
#         print 'i', i
        event.set_r(i)
        event.calc_Natural()
        int = event.get_detectableIntensity()
        xList.append(i)
        yList.append(int)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = np.array(yList)
    data['label'] = 'natural'
    
#     plot_POPS_calib([data], log=(1,1), title = 'predicted POPS response as fkt of Particle diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,event.n))
    save([data], name = 'POPSint_lamb%s_n%s'%(event.wavelength, event.n))

def calc_intensityAsFktOfRadius_test():
    """ based on: calc_intensityAsFktOfRadius_versusRuShan
        -  the aim is to get a pattern at 90 degrees and almost no area of the mirror
    """
    event = Mie(silent = True)
    # event.set_r(1.00)
    event.set_wavelength(.405)
    event.set_n(1.52)
    event.set_nang(1000)
    
    
    event.POPSdimensions['mirror diameter (mm)'] = 1
    event.POPSdimensions['mirror(top)-jet distance (mm)'] = 1000
    event.POPSdimensions['angle: jet-mirrorNormal (rad)'] = np.pi/4.
    
    xList = []
    yList = []
    for i in np.linspace(0.05,1.5,1000):
#         print 'i', i
        event.set_r(i)
        event.calc_Natural_hagen()
        int = event.get_detectableIntensity()
        xList.append(i)
        yList.append(int)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = np.array(yList)
    data['label'] = 'natural'
    
    plot_POPS_calib([data], log=(1,1), title = 'predicted POPS response as fkt of Particle diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,event.n))
#    save([data], name = 'POPSint_lamb%s_n%s'%(event.wavelength, event.n))
    
def calc_intensityAsFktOfRadius_Pops_1():
    """ creates calibration plot for the old pops"""
    event = Mie(silent = True)
    # event.set_r(1.00)
    event.set_nang(1000)
    
    
    mDiameter = event.POPSdimentions['mirror diameter (mm)'] = 25.
    rMirror = 20. # radius of the spherical mirror
    event.POPSdimentions['mirror(top)-jet distance (mm)'] = rMirror - tools.segment_hight(rMirror, mDiameter)
    event.POPSdimentions['angle: jet-mirrorNormal (rad)'] = np.pi/2.
    
    xList = []
    yList = []
    for i in np.linspace(0.05,1.5,50):
#         print 'i', i
        event.set_r(i)
        event.calc_Natural()
        int = event.get_detectableIntensity()
        xList.append(i)
        yList.append(int)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = np.array(yList)
    data['label'] = 'natural'
    
    plot_POPS_calib([data], log=(1,1), title = 'predicted POPS response as fkt of Particle diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))
#     save([data], name = 'POPSint_lamb%s_n%s'%(event.wavelength, round(event.n,2)))

def calc_intensityAsFktOfRadius_sensitivityOnRefIdx():
    """ creates calibration plot for the old pops"""
    event = Mie(silent = True, material = False)
    # event.set_r(1.00)
    event.set_nang(1000)
    
    
    event.POPSdimentions['mirror diameter (mm)'] = 25.4
    event.POPSdimentions['mirror(top)-jet distance (mm)'] = 10
    event.POPSdimentions['angle: jet-mirrorNormal (rad)'] = np.pi/2.
    
    dataList = []
    refIdxArray = np.linspace(1.59,1.63,5)
    for n in refIdxArray:
    
        xList = []
        yList = []
        event.set_n(n)
        
        for i in np.linspace(0.05,1.5,1000):
        #         print 'i', i
            event.set_r(i)
            event.calc_Natural()
            int = event.get_detectableIntensity()
            xList.append(i)
            yList.append(int)
    
        data={}
        data['x'] = 2 * np.array(xList)
        data['y'] = np.array(yList)
        data['label'] = str(round(n,4))
        dataList.append(data)
    
    plot_POPS_calib(dataList, log=(1,1), title = 'predicted POPS response as fkt of Particle diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))
    save(dataList, name = 'POPSint_lamb%s_n%s'%(event.wavelength, round(event.n,3)))

def gaussian_broadening():
    """ this is for PSLs"""
    wavelengthArray = np.linspace(.390,.420,50)
    noOfPts = np.linspace(0.05,1.5,5000)
    
    event = Mie(silent = True, design = 'POPS 1')
    event.set_nang(1000)

    dataList = []
    sumArray = np.zeros(len(noOfPts))
    for l in wavelengthArray:
    
        xList = []
        yList = []
        event.set_wavelength(l)
        
        for i in noOfPts:
        #         print 'i', i
            event.set_r(i)
            event.calc_Natural_hagen()
            int =  event.get_detectableIntensity()
            xList.append(i)
            yList.append(int)
            

    
        data={}
        data['x'] = 2 * np.array(xList)
        data['y'] = np.array(yList)
        data['label'] = str(round(l,5))+ '_' + str(round(event.n,4))
        dataList.append(data)
        
        gaussianScaling = tools.gauss_function(l, .405, 0.005) / 100.
        print(gaussianScaling)
        sumArray += gaussianScaling * data['y']
        
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray
    data['label'] = 'gaus. broadened'
    dataList.append(data)
    
    save(dataList, name = 'gaussianBroadening_')
    plot_POPS_calib(dataList, log=(1,1), title = 'predicted POPS response as fkt of Particle diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))

def dioctyl_sebacate():
    """with gaussian broadening
        - this is for Pops 1!!"""
    wavelengthArray = np.linspace(.390,.420,50)
    noOfPts = np.linspace(0.05,1.5,5000) #radius range
    
    event = Mie(silent = True, design = 'POPS 1', indexOfRef = 1.455) # Ref. for indexOfRef.  Patterson 2004
    event.set_nang(1000)

    dataList = []
    sumArray = np.zeros(len(noOfPts))
    for l in wavelengthArray:
    
        xList = []
        yList = []
        event.set_wavelength(l)
        
        for i in noOfPts:
        #         print 'i', i
            event.set_r(i)
            event.calc_Natural_hagen()
            int =  event.get_detectableIntensity()
            xList.append(i)
            yList.append(int)
            

    
        data={}
        data['x'] = 2 * np.array(xList) #this way we have a diameter range
        data['y'] = np.array(yList)
        data['label'] = str(round(l,5))+ '_' + str(round(event.n,4))
        dataList.append(data)
        
        gaussianScaling = tools.gauss_function(l, .405, 0.005) / 100.
        sumArray += gaussianScaling * data['y']
        
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray
    data['label'] = 'gaus_broadened'
    dataList.append(data)
    
    save(dataList, name = 'gaussianBroadening_')
    plot_POPS_calib(dataList, log=(1,1), title = 'predicted POPS response as fkt of Dioctyl Sebacate diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))
        

def dioctyl_sebacate_netrual_v_paraPerp():
    """ - calculates not only natural, but also parallel and perpendicular
        - with gaussian broadening
        - this is for Pops 1!!"""
    wavelengthArray = np.linspace(.390,.420,50)
    noOfPts = np.logspace(np.log10(0.05),np.log10(1.5),500) #radius range
    
    event = Mie(silent = True, design = 'POPS 1', indexOfRef = 1.455, diameter = 'dynamic') # Ref. for indexOfRef.  Patterson 2004
    event.set_nang(1000)

    dataList = []
    sumArray_natural = np.zeros(len(noOfPts))
    sumArray_parallel = np.zeros(len(noOfPts))
    sumArray_perpendicular = np.zeros(len(noOfPts))
    
    for l in wavelengthArray:
    
        xList = []
        yList_natural = []
        yList_parallel = []
        yList_perpendicular = []
        event.set_wavelength(l)
        
        for i in noOfPts:
        #         print 'i', i
            event.set_r(i)
#             event.calc_Natural_hagen()
            xList.append(i)
            
            int =  event.get_detectableIntensity("natural")
            yList_natural.append(int)
            int =  event.get_detectableIntensity("parallel")
            yList_parallel.append(int)
            int =  event.get_detectableIntensity("perpendicular")
            yList_perpendicular.append(int)
            

    
#         data={}
#         data['x'] = 2 * np.array(xList) #this way we have a diameter range
#         data['y'] = np.array(yList)
#         data['label'] = str(round(l,5))+ '_' + str(round(event.n,4))
#         dataList.append(data)
        
        gaussianScaling = tools.gauss_function(l, .405, 0.005) / 100.
        sumArray_natural += gaussianScaling * np.array(yList_natural)
        sumArray_parallel += gaussianScaling * np.array(yList_parallel)
        sumArray_perpendicular += gaussianScaling * np.array(yList_perpendicular)
        
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_natural
    data['label'] = 'natural'
    dataList.append(data)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_parallel
    data['label'] = 'parallel'
    dataList.append(data)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_perpendicular
    data['label'] = 'perpendicular'
    dataList.append(data)
    
    save(dataList, name = '')
    plot_POPS_calib(dataList, log=(1,1), title = 'predicted POPS response as fkt of Dioctyl Sebacate diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))
        
def dioctyl_sebacate_netrual_v_paraPerp_exactAngleDependence():
    """ - calculates not only natural, but also parallel and perpendicular
        - with gaussian broadening
        - this is for Pops 1!!"""
    wavelengthArray = np.linspace(.400,.410,15)
    noOfPts = np.logspace(np.log10(0.05),np.log10(1.5),30) #radius range 
    
    event = Mie(silent = True, design = 'POPS 2', indexOfRef = 1.455, diameter = 'dynamic') # Ref. for indexOfRef.  Patterson 2004
    event.set_nang(20)

    dataList = []
    sumArray_natural = np.zeros(len(noOfPts))
    sumArray_parallel = np.zeros(len(noOfPts))
    sumArray_perpendicular = np.zeros(len(noOfPts))
    
    for e,l in enumerate(wavelengthArray):
        print('start wavelength %s of %s' % (e, len(wavelengthArray)))
        xList = []
        yList_natural = []
        yList_parallel = []
        yList_perpendicular = []
        event.set_wavelength(l)
        
        for i in noOfPts:
#             print 'start a point'
        #         print 'i', i
            event.set_r(i)
#             event.calc_Natural_hagen()
            xList.append(i)
            int =  event.get_detectableIntensity("natural")
            yList_natural.append(int)
            int =  event.get_detectableIntensity("parallel")
            yList_parallel.append(int)
            int =  event.get_detectableIntensity("perpendicular")
            yList_perpendicular.append(int)
            

    
#         data={}
#         data['x'] = 2 * np.array(xList) #this way we have a diameter range
#         data['y'] = np.array(yList)
#         data['label'] = str(round(l,5))+ '_' + str(round(event.n,4))
#         dataList.append(data)
        
        gaussianScaling = tools.gauss_function(l, .405, 0.005) / 100.
        sumArray_natural += gaussianScaling * np.array(yList_natural)
        sumArray_parallel += gaussianScaling * np.array(yList_parallel)
        sumArray_perpendicular += gaussianScaling * np.array(yList_perpendicular)
        
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_natural
    data['label'] = 'natural'
    dataList.append(data)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_parallel
    data['label'] = 'parallel'
    dataList.append(data)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_perpendicular
    data['label'] = 'perpendicular'
    dataList.append(data)
    
    save(dataList, name = '')
    plot_POPS_calib(dataList, log=(1,1), title = 'predicted POPS response as fkt of Dioctyl Sebacate diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))
            
def dioctyl_sebacate_variousIndices():
    """ - calculates not only natural, but also parallel and perpendicular
        - with gaussian broadening
        - this is for Pops 1!!"""
    wavelengthArray = np.linspace(.390,.420,30)
    noOfPts = np.logspace(np.log10(0.05),np.log10(1.5),200) #radius range
    
    event = Mie(silent = True, design = 'POPS 2', indexOfRef = 1.455, diameter = 'dynamic') # Ref. for indexOfRef.  Patterson 2004
    event.set_nang(100)

    dataList = []
    sumArray_natural = np.zeros(len(noOfPts))
    sumArray_parallel = np.zeros(len(noOfPts))
    sumArray_perpendicular = np.zeros(len(noOfPts))
    
    for e,l in enumerate(wavelengthArray):
        print('start wavelength %s of %s' % (e, len(wavelengthArray)))
        xList = []
        yList_natural = []
        yList_parallel = []
        yList_perpendicular = []
        event.set_wavelength(l)
        
        for i in noOfPts:
#             print 'start a point'
        #         print 'i', i
            event.set_r(i)
#             event.calc_Natural_hagen()
            xList.append(i)
            int =  event.get_detectableIntensity("natural")
            yList_natural.append(int)
            int =  event.get_detectableIntensity("parallel")
            yList_parallel.append(int)
            int =  event.get_detectableIntensity("perpendicular")
            yList_perpendicular.append(int)
            

    
#         data={}
#         data['x'] = 2 * np.array(xList) #this way we have a diameter range
#         data['y'] = np.array(yList)
#         data['label'] = str(round(l,5))+ '_' + str(round(event.n,4))
#         dataList.append(data)
        
        gaussianScaling = tools.gauss_function(l, .405, 0.005) / 100.
        sumArray_natural += gaussianScaling * np.array(yList_natural)
        sumArray_parallel += gaussianScaling * np.array(yList_parallel)
        sumArray_perpendicular += gaussianScaling * np.array(yList_perpendicular)
        
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_natural
    data['label'] = 'natural'
    dataList.append(data)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_parallel
    data['label'] = 'parallel'
    dataList.append(data)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_perpendicular
    data['label'] = 'perpendicular'
    dataList.append(data)
    
    save(dataList, name = '')
    plot_POPS_calib(dataList, log=(1,1), title = 'predicted POPS response as fkt of Dioctyl Sebacate diameter ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))

def different_MirrorDistance():
    """ - calculates not only natural, but also parallel and perpendicular
        - with gaussian broadening
        - this is for Pops 1!!"""
    wavelengthArray = np.linspace(.400,.410,15)
    noOfPts = np.logspace(np.log10(0.05),np.log10(1.5),500) #radius range 
    
    event = Mie(silent = True, design = 'POPS 2', indexOfRef = 1.455, diameter = 'dynamic') # Ref. for indexOfRef.  Patterson 2004
    event.set_nang(100) # number of scatternig angles
    event.POPSdimensions['mirror(top)-jet distance (mm)'] = 7.68#+2.159
    
    dataList = []
    sumArray_natural = np.zeros(len(noOfPts))
    sumArray_parallel = np.zeros(len(noOfPts))
    sumArray_perpendicular = np.zeros(len(noOfPts))
    
    for e,l in enumerate(wavelengthArray):
        print('start wavelength %s of %s' % (e, len(wavelengthArray)))
        xList = []
        yList_natural = []
        yList_parallel = []
        yList_perpendicular = []
        event.set_wavelength(l)
        
        for i in noOfPts:
#             print 'start a point'
        #         print 'i', i
            event.set_r(i)
#             event.calc_Natural_hagen()
            xList.append(i)
            int =  event.get_detectableIntensity("natural")
            yList_natural.append(int)
            int =  event.get_detectableIntensity("parallel")
            yList_parallel.append(int)
            int =  event.get_detectableIntensity("perpendicular")
            yList_perpendicular.append(int)
            

    
#         data={}
#         data['x'] = 2 * np.array(xList) #this way we have a diameter range
#         data['y'] = np.array(yList)
#         data['label'] = str(round(l,5))+ '_' + str(round(event.n,4))
#         dataList.append(data)
        
        gaussianScaling = tools.gauss_function(l, .405, 0.005) / 100.
        sumArray_natural += gaussianScaling * np.array(yList_natural)
        sumArray_parallel += gaussianScaling * np.array(yList_parallel)
        sumArray_perpendicular += gaussianScaling * np.array(yList_perpendicular)
        
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_natural
    data['label'] = 'natural'
    dataList.append(data)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_parallel
    data['label'] = 'parallel'
    dataList.append(data)
    
    data={}
    data['x'] = 2 * np.array(xList)
    data['y'] = sumArray_perpendicular
    data['label'] = 'perpendicular'
    dataList.append(data)
    
    save(dataList, name = '14_04_18_differentMirrorDistances_%s'%event.POPSdimensions['mirror(top)-jet distance (mm)'])
    plot_POPS_calib(dataList, log=(1,1), title = 'Different mirror distances')# ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))
    return dataList


def wavelengthDependence(nm):
    dRange = np.logspace(np.log10(0.05),np.log10(1.5),200) #radius range 
    
    event = Mie(silent = True, design = 'POPS 2', indexOfRef = 1.455, diameter = 'dynamic') # Ref. for indexOfRef.  Patterson 2004
    event.set_nang(100) # number of scatternig angles
    event.POPSdimensions['mirror(top)-jet distance (mm)'] = 7.68#+2.159
    
    dataList = []
#    sumArray_natural = np.zeros(len(noOfPts))
#    sumArray_parallel = np.zeros(len(noOfPts))
#    sumArray_perpendicular = np.zeros(len(noOfPts))
    
#    for e,l in enumerate(wavelengthArray):
#        print 'start wavelength %s of %s'%(e, len(wavelengthArray))
#        xList = []
#        yList_natural = []
#        yList_parallel = []
#        yList_perpendicular = []
    #nm = 405
    event.set_wavelength(nm)
    perpInt = []
    for i in dRange:
#             print 'start a point'
    #         print 'i', i
        event.set_r(i)
#             event.calc_Natural_hagen()

        int =  event.get_detectableIntensity("perpendicular")
        perpInt.append(int)
            
        
    data={}
    data['x'] = 2 * np.array(dRange)
    data['y'] = perpInt
    data['label'] = '%s nm'%nm
    dataList.append(data)

    
    save(dataList, name = '14_05_14_nmDependence_IOR1.455')
    plot_POPS_calib(dataList, log=(1,1), title = 'Different mirror distances')# ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))
    return dataList
    
def DOSversusNHSO(IOR):
    """dioctyl sebacate(1.455) versus Ammonium sulfate (1.40)"""
    dRange = np.logspace(np.log10(0.05),np.log10(1.5),200) #radius range 
    
    event = Mie(silent = True, design = 'POPS 2', indexOfRef = IOR, diameter = 'dynamic') # Ref. for indexOfRef.  Patterson 2004
    event.set_nang(100) # number of scatternig angles
    event.POPSdimensions['mirror(top)-jet distance (mm)'] = 7.68#+2.159
    
    dataList = []
#    sumArray_natural = np.zeros(len(noOfPts))
#    sumArray_parallel = np.zeros(len(noOfPts))
#    sumArray_perpendicular = np.zeros(len(noOfPts))
    
#    for e,l in enumerate(wavelengthArray):
#        print 'start wavelength %s of %s'%(e, len(wavelengthArray))
#        xList = []
#        yList_natural = []
#        yList_parallel = []
#        yList_perpendicular = []
    nm = .405
    event.set_wavelength(nm)
    perpInt = []
    for i in dRange:
#             print 'start a point'
    #         print 'i', i
        event.set_r(i)
#             event.calc_Natural_hagen()

        int =  event.get_detectableIntensity("perpendicular")
        perpInt.append(int)
            
        
    data={}
    data['x'] = 2 * np.array(dRange)
    data['y'] = perpInt
    data['label'] = '%s ior'%IOR
    dataList.append(data)

    
    save(dataList, name = '14_05_16_IOR_dipendenc_405nm')
    plot_POPS_calib(dataList, log=(1,1), title = 'Different mirror distances')# ($\lambda=$%s nm; $n=$%s)'%(event.wavelength,round(event.n,2)))
    return dataList
##################################################################################################################################################################
##################################################################################################################################################################
if __name__ == "__main__":
    print('los gehts')
#    test_comparison_to_internet()
#     test_calc_nintyOnly()
#     default()
#     bhmie_hagen_test()
#     calc_intensityAsFktOfRadius_sensitivityOnRefIdx()
#     calc_scatteringPattern()
#     calc_nintyOnly()
#    gaussian_broadening()
#    dioctyl_sebacate()
#     dioctyl_sebacate_netrual_v_paraPerp()
#    calc_intensityAsFktOfRadius_test()
#    dioctyl_sebacate_netrual_v_paraPerp_exactAngleDependence()
#     dioctyl_sebacate_variousIndices()



