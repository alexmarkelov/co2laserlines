# Programm make array of points associated with envelope
# curve of amplification line of CO2 molecular

import math
import numpy as np

PLANK = 6.62607015e-34
LIGHTVEL = 299792458
KBOLTZMANN = 1.380649e-23
PI = 3.14159265358979323846264
ROTCONST = 38.7141392007267
RGAS = 8.31431
NAVOGADRO = 6.022140*10**23


class LineWidthConst:
    

def linewidthconst(branch, wavelength, jnum):
    width = None
    if branch == 'P' and wavelength == 10:
        width = (-0.0733*jnum + 9.1664)/2
    elif branch == 'R' and wavelength == 10:
        width = (-0.0771*jnum + 9.0874)/2
    elif branch == 'P' and wavelength == 9:
        width = (-0.0563*jnum + 8.7477)/2
    elif branch == 'R' and wavelength == 9:
        width = (-0.0479*jnum + 8.1615)/2
    else:
        width = None
    return width

def getlinewidth(lwconst, pressure, temp, CO2, N2, He, Xe):
    lwidth = lwconst*10e6/PI*pressure*(CO2 + 0.73*N2 + 0.64*(He + Xe))*math.sqrt(300/temp)
    return lwidth

def lineform(linewidth, wl, frequency):
    return (linewidth/2*PI)*((frequency - wl)**2+(linewidth/2)**2)**-1

def ngasmixture(CO2, N2, He, Xe):
    return (1 - Xe)*(CO2*44 + N2*28 + He*4) + Xe*131.29

def numberofCO2mol(weight, heght, length, pressure, temp, CO2, He, Ne ):
    V = weight*heght*length
    molarmass = ngasmixture(CO2, N2, He, Ne)
    ng = (pressure*molarmass*V)/(RGAS*temp)
    numberCO2mol = (CO2*44/molarmass)*ng*NAVOGADRO
    return numberCO2mol

def coefB(branch, wl, jnum):
    if branch == 'P' and wl == 10:
        B = -0.0009*jnum + 0.1902
    elif branch == 'R' and wl == 10:
        B = -0.001*jnum + 0.2035
    elif branch == 'P' and wl == 9:
        B = 0.0011*jnum + 0.1542
    elif branch == 'R' and wl == 9:
        B = 0.0014*jnum + 0.1567
    else:
        B = None
    return B

def energyF(self, j):
    return ROTCONST*j*(j + 1)

    
class CO2AmplificationCurve:

    def __init__(self, jnum, branch, central_wavelength, pressure,
                 temperature, CO2part, N2part, Hepart, Xepart, df, Nmax,
                 laser_width, laser_heigth, laser_length, pump_coeff):
        self.jnum  = jnum
        self.branch = branch
        self.wl = central_wavelength
        self.pr = pressure
        self.temp = temperature
        self.CO2 = CO2part
        self.N2 = N2part
        self.He = Hepart
        self.Xe = Xepart
        self.df = df
        self.Nmax = int(Nmax)
        self.lasw = laser_width
        self.lash = laser_heigth
        self.lasl = laser_length
        self.pump = pump_coeff
        self.ampcurve = np.zeros((self.Nmax, 2))
        self._lwconst = linewidthconst(self.branch, self.wl, self.jnum)
        self._linewidth = getlinewidth(self._lwconst, self.pr, self.temp, self.CO2, self.N2, self.He, self.Xe)
        self._numCO2mol = numberofCO2mol()
        getlineform()

    

    
        

    





    def constAmp(self):
        B = coefB()
        F = energyF(self.jnum)
        E = PLANK*LIGHTVEL/(KBOLTZMANN*self.temp)
        С = ((LIGHTVEL/self.wl)**2)*PLANK*LIGHTVEL*ROTCONST/(4*PI*KBOLTZMANN*self.temp)
        Ndown = self._numCO2mol*math.exp(-PLANK*LIGHTVEL*100*1400/(KBOLTZMANN*self.temp))
        Nup = self._numCO2mol*self.pump - Ndown
        if self.branch == 'P' and self.wl == 10:
            D = C*(2*self.jnum - 1)*(coefB())*(Nup*math.exp(-energyF(self.jnum - 1)*E)
                                           - Ndown*math.exp(-energyF(self.jnum)*E))
            return D
        elif self.branch == 'R' and self.wl == 10:
            D = C*(2*self.jnum - 1)*(coefB())*(Nup*math.exp(-energyF(self.jnum - 1)*E)
                                           - Ndown*math.exp(-energyF(self.jnum)*E))
            return D
        elif self.branch == 'P' and self.wl == 9:
            D = C*(2*self.jnum + 3)*(coefB())*(Nup*math.exp(-energyF(self.jnum + 1)*E)
                                           - Ndown*math.exp(-energyF(self.jnum)*E))
            return D
        elif self.branch == 'P' and self.wl == 9:
            D = C*(2*self.jnum + 3)*(coefB())*(Nup*math.exp(-energyF(self.jnum + 1)*E)
                                           - Ndown*math.exp(-energyF(self.jnum)*E))
            return D
        else:
            return None
        
    def getlineform(self):
        for index in range(self.Nmax):
            freq = self.wl - self.df*(int(self.Nmax/2)) + self.df*index
            self.ampcurve[index, 0] = freq
            self.ampcurve[index, 1] = constAmp()*lineform(frequency)
        pass
if __name__ == "__main__":
    A = CO2AmplificationCurve(20, 'P', 10, 160*133.32, 400, 1, 1, 5, 0.05,
                          1*10e5, 10e4, 2*10e-3, 2*10e-3, 150*10e-3, 0.25)
    

    
