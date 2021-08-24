# Programm make array of points associated with envelope
# curve of amplification line of CO2 molecular

import math
import numpy as np
import CO2lines_wavelength as wl

PLANK = 6.62607015e-34
LIGHT_VELOCITY = wl.LIGTH_VELOCITY
KBOLTZMANN = 1.380649e-23
PI = 3.14159265358979323846264
ROTCONST = 38.7141392007267
RGAS = 8.31431
NAVOGADRO = 6.022140*10**23


line_width_const = {
    '10P': 'line_width_const_10P',
    '10R': 'line_width_const_10R',
    '9P': 'line_width_const_9P',
    '9R': 'line_width_const_9R'
    }
coef_B_const = {
    '10P': 'coef_B_const_10P',
    '10R': 'coef_B_const_10R',
    '9P': 'coef_B_const_9P',
    '9R': 'coef_B_const_9R'
    }



def line_width_const_10P(jnum: int) -> float:
    width = (-0.0733*jnum + 9.1664)/2
    return width


def line_width_const_10R(jnum: int) -> float:
    width = (-0.0771*jnum + 9.0874)/2
    return width


def line_width_const_9P(jnum: int) -> float:
    width = (-0.0563*jnum + 8.7477)/2
    return width


def line_width_const_9R(jnum: int) -> float:
    width = (-0.0479*jnum + 8.1615)/2
    return width


def coef_B_const_10P(jnum: int) -> float:
    B = -0.0009*jnum + 0.1902
    return B


def coef_B_const_10R(jnum: int) -> float:
    B = -0.001*jnum + 0.2035
    return B


def coef_B_const_9P(jnum: int) -> float:
    B = 0.0011*jnum + 0.1542
    return B


def coef_B_const_9R(jnum: int) -> float:
    B = 0.0014*jnum + 0.1567
    return B


def get_line_coef_B(line: str) -> float:
    center_branch = wl.get_line_center_branch(line)
    jnum = wl.get_line_jnum(line)
    const_B_calc = eval(coef_B_const[center_branch])
    B = const_B_calc(jnum)
    return B


assert(get_line_coef_B('10_P_20'))


def get_line_width_const(line: str) -> float:
    center_branch = wl.get_line_center_branch(line)
    jnum = wl.get_line_jnum(line)
    line_width_calc = eval(line_width_const[center_branch])
    width = line_width_calc(jnum)
    return width


assert(get_line_width_const('10_P_20'))


class GasMixture:
    def __init__(self, *args, pressure=100, temp=300,
                 CO2=1, N2=1, He=5, Xe=0.05):
        self.pressure = pressure
        self.temp = temp
        self.CO2 = CO2/(CO2 + N2 + He +Xe*(CO2 + N2 + He))
        self.N2 = N2/(CO2 + N2 + He +Xe*(CO2 + N2 + He))
        self.He = He/(CO2 + N2 + He +Xe*(CO2 + N2 + He))
        self.Xe = Xe/(CO2 + N2 + He +Xe*(CO2 + N2 + He))
        self.n_gas = self.get_n_gas_mixture()

    def get_n_gas_mixture(self) -> float:
        n_gas = (1 - self.Xe)*(self.CO2*44 + self.N2*28 +
                               self.He*4) + self.Xe*131.29
        return n_gas


def get_line_width(line: str, gas_mixture: GasMixture = GasMixture() ):
    line_width_const = get_line_width_const(line)
    line_width = line_width_const*10e6*gas_mixture.pressure*(
        gas_mixture.CO2 + 0.73*gas_mixture.N2 + 0.64*(
            gas_mixture.He + gas_mixture.Xe))*math.sqrt(300/gas_mixture.temp)
    return line_width


def lineform(linewidth, wl, frequency):
    return (linewidth/(2*PI))/((frequency - wl)**2+(linewidth/2)**2)

gas_mixture_1 = GasMixture()

class GasVolume:
    def __init__(self, *args, width=0.3, height=0.3, length=13.6, diameter=None):
        if diameter:
            self.length = length
            self.diameter = diametr
            self.volume = PI*diameter**2/4*length
        else:
            self.length = length
            self.width = width
            self.height = height
            self.volume = length*width*height


def number_CO2_mol(gas_volume=GasVolume(), gas_mixture=GasMixture()):
    mol_gas = (gas_mixture.pressure*gas_mixture.n_gas*gas_volume.volume)/(
        RGAS*temp)
    number_CO2_mol = (CO2*44/gas_mixture.n_gas)*mol_gas*NAVOGADRO
    return number_CO2_mol


def energy_F(line: str) -> float:
    jnum = wl.get_line_jnum(line)
    return ROTCONST*jnum*(jnum + 1)


def member_A(line: str, gas_mixture: GasMixture):
    wl = 
    ((LIGHTVEL/wl)**2)*PLANK*LIGHTVEL*ROTCONST/(4*PI*KBOLTZMANN*temp)


class CO2AmplificationCurve:
    def __init__(self, line, *args, gas_mixture=GasMixture(), df=1e5, Nmax=10000,
                 gas_volume=GasVolume, pump_coeff=0.5):
        self.line  = line
        self.gas_mixture = gas_mixture
        self.gas_volume = gas_volume
        self.pump_coeff = pump_coeff
        self.df = df
        self.Nmax = Nmax
        self.amp_curve = np.zeros((self.Nmax, 2))
        self.line_ = get_line_width(self.line, self.gas_mixture)
        self.numCO2mol = number_CO2_mol(self.gas_volume, self.gas_mixture)
    
    def const_ampl(self):
        B = get_line_coef_B(self.line)
        F = energyF(self.jnum)
        E = PLANK*LIGHTVEL/(KBOLTZMANN*self.temp)
        С = ((LIGHTVEL/self.wl)**2)*PLANK*LIGHTVEL*ROTCONST/(4*PI*KBOLTZMANN*self.temp)
        Ndown = self.numCO2mol*math.exp(-PLANK*LIGHTVEL*100*1400/(KBOLTZMANN*self.temp))
        Nup = self.numCO2mol*self.pump - Ndown
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
    

    
