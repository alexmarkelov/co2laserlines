# Programm make array of points associated with envelope
# curve of amplification line of CO2 molecular

import math
import numpy as np
import CO2lines_wavelength as wl

PLANK = 6.62607015e-34
LIGHT_VELOCITY = wl.LIGHT_VELOCITY
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
coef_number_degeneration = {
    '10': 'number_degeneration_10mkm',
    '9': 'number_degeneration_9mkm'
    }

up_line_branch_const ={
    'P': 'self.line.jnum-1',
    'R': 'self.line.jnum+1'
    }

def up_line_branch(line: wl.Line) -> int:
    return up_line_branch_const[line.branch]

def number_degeneration_10mkm(jnum: int):
    return 2*jnum-1

def number_degeneration_9mkm(jnum: int):
    return 2*jnum+3

def number_degeneration_calc(line: wl.Line):
    return eval(coef_number_degeneration[str(line.center)])(line.jnum)

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


def get_line_coef_B(line: wl.Line) -> float:
    center_branch = str(line.center) + line.branch
    jnum = line.jnum
    const_B_calc = eval(coef_B_const[center_branch])
    B = const_B_calc(jnum)
    return B

line10P20 = wl.Line(10, 'P', 20)
assert(get_line_coef_B(line10P20))


def get_line_width_const(line: wl.Line) -> float:
    center_branch = str(line.center) + line.branch
    jnum = line.jnum
    line_width_calc = eval(line_width_const[center_branch])
    width = line_width_calc(jnum)
    return width


assert(get_line_width_const(line10P20))


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


def lineform(frequency: float, line_center_frequency: float,
             half_line_width: float):
    return (half_line_width/PI)/((frequency - line_center_frequency)**2+
                                  half_line_width**2)

gas_mixture_1 = GasMixture()

class GasVolume:
    def __init__(self, *args, width=0.3, height=0.3, length=13.6, diameter=None):
        if diameter:
            self.length = length
            self.diameter = diametr
            self.volume = PI*diameter**2/4*length/1000000
        else:
            self.length = length
            self.width = width
            self.height = height
            self.volume = length*width*height/1000000

gas_volume1  = GasVolume()
#print(gas_volume1.volume)


def number_CO2_mol(gas_volume=GasVolume(), gas_mixture=GasMixture()):
    mol_gas = (gas_mixture.pressure*gas_mixture.n_gas*gas_volume.volume)/(
        RGAS*gas_mixture.temp)
    number_CO2_mol = (gas_mixture.CO2*44/gas_mixture.n_gas)*mol_gas*NAVOGADRO
    return number_CO2_mol


def get_energy_F(jnum) -> float:
    return ROTCONST*jnum*(jnum + 1)


def member_A(line: wl.Line, gas_mixture: GasMixture):
    wl = line.wavelength
    value =(LIGHT_VELOCITY/wl)**2*PLANK*LIGHT_VELOCITY*ROTCONST/(
        4*PI*KBOLTZMANN*gas_mixture.temp)
    return value

def member_B(gas_mixture: GasMixture):
    return PLANK*LIGHT_VELOCITY/(KBOLTZMANN*gas_mixture.temp)


def number_up_down_mol(number_CO2_mol: float, gas_mixture: GasMixture,
                       pump_coeff: float=0.5):
    Ndown = number_CO2_mol*math.exp(-PLANK*LIGHT_VELOCITY*100*1400/(
        KBOLTZMANN*gas_mixture.temp))
    Nup = number_CO2_mol*pump_coeff - Ndown
    return {'Nup': Nup, 'Ndown': Ndown}

class CO2AmplificationCurve:
    def __init__(self, line, *args, gas_mixture=GasMixture(), df=1e5,
                 Nmax=10000, gas_volume=GasVolume(), pump_coeff=0.5):
        self.line  = line
        self.gas_mixture = gas_mixture
        self.gas_volume = gas_volume
        self.pump_coeff = pump_coeff
        self.df = df
        self.Nmax = Nmax
        self.amp_curve_freq = list()
        self.amp_curve_ampl = list()
        self.amp_curve = np.zeros((self.Nmax*2, 2))
        self.half_line_width = get_line_width(self.line, self.gas_mixture)/2
        self.numCO2mol = number_CO2_mol(self.gas_volume, self.gas_mixture)
        self.amp_coef = self.const_ampl()
        self.get_line_form()
    
    def const_ampl(self):
        A = member_A(self.line, self.gas_mixture)
        #print(A)
        B = member_B(self.gas_mixture)
        #print(B)
        coef_B = get_line_coef_B(self.line)
        #print(coef_B)
        up_down = number_up_down_mol(self.numCO2mol,
                                            self.gas_mixture, self.pump_coeff)
        number_degeneration = number_degeneration_calc(self.line)
        #print(up_down, number_degeneration)
        const_ampl = A*coef_B*number_degeneration_calc(self.line)*(
            up_down['Nup']*math.exp(-get_energy_F(eval(
                up_line_branch(self.line)))*B)- up_down['Ndown']*math.exp(
                    -get_energy_F(self.line.jnum)*B))
        return const_ampl
        
        
    def get_line_form(self):
        for index in range(2*self.Nmax + 1):
            freq = self.line.wavelength + self.df*self.Nmax - self.df*index
            self.amp_curve_freq.append(freq)
            self.amp_curve_ampl.append(self.amp_coef*lineform(freq,
                            self.line.wavelength, self.half_line_width))
        pass

if __name__ == "__main__":
    line10P20 = wl.Line(10,'P', 20, unit='hz' )
    gas_mixture1 = GasMixture(temp=400)
    gas_volume1 = GasVolume()
    branch_10P = wl.Branch(10, 'P', 2, 58, unit='hz')
    for line in branch_10P:
        amplification = CO2AmplificationCurve(line, gas_mixture = gas_mixture1, gas_volume = gas_volume1)
        print(line.name, ' - ',amplification.amp_coef)


    
