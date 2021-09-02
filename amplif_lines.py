'''
Module intended to calculate a curve of amplification line of CO2 molecular.
Four branches (R, P on 10,9 mkm) are available.

amplif_lines includes 3 classes:

GasMixture is a class containing information about composition, temperature
and pressure of gas mixture. Example of using:
GasMixture(pressure=100, temp=350, CO2=1, N2=1, He=5, Xe=0.05)
where pressure in torr, temp in grad K. CO2=1, N2=1, He=5, Xe=0.05 means
gas mixture composition CO2:N2:He = 1:1:5 volume fractions plus 5% volume
fraction of Xe.

GasVolume is a class containin information about dimensions in cm of volume
with gas. Example of using:
GasMixture(width=1, height=1, length=15) for rectangular parallelepiped volume
and GasMixture(length=15, diameter=0.5) for cylindrical volume. 

AmplificationCurve is a class containing information about amplification curve
of one emission line of 12C16O2 molecule. Example of using:
AmplificationCurve(line, gas_mixture = gas_mixture_1,
gas_volume = gas_volume_1, df=1e5, Nmax=10000, gas_volume=GasVolume(),
pump_coeff=0.5),where:
line is emission line, Line class instance from wavelength module,
gas_mixture_1 is GasMixture class instance,
df is frequency step,
2*Nmax + 1 is a number of steps, for computing CO2 line amplification curve.
''' 


import math
import wavelength as wl

PLANK = 6.62607015e-34
LIGHT_VELOCITY = wl.LIGHT_VELOCITY
KBOLTZMANN = 1.380649e-23
PI = math.pi
ROTCONST = 38.7141392007267
RGAS = 8.31431
NAVOGADRO = 6.022140*10**23

#The dictionaries bellow are needed for calculation amlification curve 
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
    '''Return jnum of upper level'''
    return up_line_branch_const[line.branch]


def number_degeneration_10mkm(jnum: int):
    '''Return the degeneration of lower level for 10 mkm branch'''
    return 2*jnum-1


def number_degeneration_9mkm(jnum: int):
    '''Return the degeneration of lower level for 9 mkm branch'''
    return 2*jnum+3

def number_degeneration_calc(line: wl.Line):
    '''Return the degeneration for the line'''
    return eval(coef_number_degeneration[str(line.center)])(line.jnum)

def line_width_const_10P(jnum: int) -> float:
    '''Return width for 10P branch with J = jnum'''
    width = (-0.0733*jnum + 9.1664)/2
    return width


def line_width_const_10R(jnum: int) -> float:
    '''Return width for 10R branch with J = jnum'''
    width = (-0.0771*jnum + 9.0874)/2
    return width


def line_width_const_9P(jnum: int) -> float:
    '''Return width for 9P branch with J = jnum'''
    width = (-0.0563*jnum + 8.7477)/2
    return width


def line_width_const_9R(jnum: int) -> float:
    '''Return width for 9R branch with J = jnum'''
    width = (-0.0479*jnum + 8.1615)/2
    return width


def coef_B_const_10P(jnum: int) -> float:
    '''Return Einstein's coefficient B for 10P branch with J = jnum'''
    B = -0.0009*jnum + 0.1902
    return B


def coef_B_const_10R(jnum: int) -> float:
    '''Return Einstein's coefficient B for 10R branch with J = jnum'''
    B = -0.001*jnum + 0.2035
    return B


def coef_B_const_9P(jnum: int) -> float:
    '''Return Einstein's coefficient B for 9P branch with J = jnum'''
    B = 0.0011*jnum + 0.1542
    return B


def coef_B_const_9R(jnum: int) -> float:
    '''Return Einstein's coefficient B for 9R branch with J = jnum'''
    B = 0.0014*jnum + 0.1567
    return B


def get_line_coef_B(line: wl.Line) -> float:
    '''Return Einstein's coefficient B for line'''
    center_branch = str(line.center) + line.branch
    jnum = line.jnum
    const_B_calc = eval(coef_B_const[center_branch])
    B = const_B_calc(jnum)
    return B


def get_line_width_const(line: wl.Line) -> float:
    '''Return width amplification curve for line'''
    center_branch = str(line.center) + line.branch
    jnum = line.jnum
    line_width_calc = eval(line_width_const[center_branch])
    width = line_width_calc(jnum)
    return width


class GasMixture:
    '''GasMixture is a class for description of gas mixture parameters,
       parameters: pressure[60-200 torr], temp(temperature)[150-500 grad K],
       gas mixture composition - CO2, N2, He - volume part of the gas
       component, for example 1:1:5 = CO2, N2, He;
       Xe - percentage of gas component[rel. unit]'''
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
    '''Return line width for input line and gas mixture'''
    line_width_const = get_line_width_const(line)
    line_width = line_width_const*10e6*gas_mixture.pressure*(
        gas_mixture.CO2 + 0.73*gas_mixture.N2 + 0.64*(
            gas_mixture.He + gas_mixture.Xe))*math.sqrt(300/gas_mixture.temp)
    return line_width


def lineform(frequency: float, line_center_frequency: float,
             half_line_width: float):
    '''The function determine the line form'''
    return (half_line_width/PI)/((frequency - line_center_frequency)**2+
                                  half_line_width**2)


class GasVolume:
    '''The class defining the volume of acitve gas medium of laser'''
    def __init__(self, *args, width=0.3, height=0.3, length=13.6,
                 diameter=None):
        if diameter:
            self.length = length
            self.diameter = diametr
            self.volume = PI*diameter**2/4*length/1000000
        else:
            self.length = length
            self.width = width
            self.height = height
            self.volume = length*width*height/1000000


def number_CO2_mol(gas_volume=GasVolume(), gas_mixture=GasMixture()):
    '''Return the number of CO2 molecules in gas volume
       with parmeter of gas mixture '''
    mol_gas = (gas_mixture.pressure*gas_mixture.n_gas*gas_volume.volume)/(
        RGAS*gas_mixture.temp)
    number_CO2_mol = (gas_mixture.CO2*44/gas_mixture.n_gas)*mol_gas*NAVOGADRO
    return number_CO2_mol


def get_energy_F(jnum) -> float:
    '''Return the rotational energy for jnum'''
    return ROTCONST*jnum*(jnum + 1)


def member_A(line: wl.Line, gas_mixture: GasMixture):
    '''Get parameter A for amplificarion line magnitude,
       depends on line and gas mixture'''
    wl = line.wavelength
    value =(LIGHT_VELOCITY/wl)**2*PLANK*LIGHT_VELOCITY*ROTCONST/(
        4*PI*KBOLTZMANN*gas_mixture.temp)
    return value

def member_B(gas_mixture: GasMixture):
    '''Get parameter B for amplification line magnitude,
       depends on gas mixture'''
    return PLANK*LIGHT_VELOCITY/(KBOLTZMANN*gas_mixture.temp)


def number_up_down_mol(number_CO2_mol: float, gas_mixture: GasMixture,
                       pump_coeff: float=0.5):
    '''Returns the number of CO2 molecules on upper and lower level,
       depends on number of all CO2 molecules, gas mixture and pump
       coefficient'''
    Ndown = number_CO2_mol*math.exp(-PLANK*LIGHT_VELOCITY*100*1400/(
        KBOLTZMANN*gas_mixture.temp))
    Nup = number_CO2_mol*pump_coeff - Ndown
    return {'Nup': Nup, 'Ndown': Ndown}

class AmplificationCurve:
    def __init__(self, line: wl.Line, *args, gas_mixture=GasMixture(), df=1e5,
                 Nmax=10000, gas_volume=GasVolume(), pump_coeff=0.5):
        self.__line  = line
        self.__gas_mixture = gas_mixture
        self.__gas_volume = gas_volume
        self.__pump_coeff = pump_coeff
        self.__df = df
        self.__Nmax = Nmax
        self.__amp_curve_freq = list()
        self.__amp_curve_ampl = list()
        self.__half_line_width = get_line_width(self.__line,
                                                self.__gas_mixture)/2
        self.__numCO2mol = number_CO2_mol(self.__gas_volume,
                                          self.__gas_mixture)
        self.get_const_ampl()
        self.get_line_form()

    @property
    def line(self):
        return self.__line

    @property
    def gas_mixture(self):
        return self.__gas_mixture

    @property
    def gas_volume(self):
        return self.__gas_volume

    @property
    def pump_coeff(self):
        return self.__pump_coeff

    @pump_coeff.setter
    def pump_coeff(self, new_pump_coeff):
        self.__pump_coeff = new_pump_coeff
        self.get_const_ampl()
        self.get_line_form()
        return new_pump_coeff

    @property
    def amplif(self):
        return self.__amp_curve_ampl[self.__Nmax]

    @property
    def amp_curve_freq(self):
        return self.__amp_curve_freq

    @property
    def amp_curve_ampl(self):
        return self.__amp_curve_ampl
    
    def get_const_ampl(self):
        A = member_A(self.__line, self.__gas_mixture)
        #print(A)
        B = member_B(self.__gas_mixture)
        #print(B)
        coef_B = get_line_coef_B(self.__line)
        #print(coef_B)
        up_down = number_up_down_mol(self.__numCO2mol,
                                            self.__gas_mixture, self.__pump_coeff)
        number_degeneration = number_degeneration_calc(self.__line)
        #print(up_down, number_degeneration)
        const_ampl = A*coef_B*number_degeneration_calc(self.__line)*(
            up_down['Nup']*math.exp(-get_energy_F(eval(
                up_line_branch(self.__line)))*B)- up_down['Ndown']*math.exp(
                    -get_energy_F(self.__line.jnum)*B))
        self.__amp_coeff = const_ampl
        return const_ampl
        
        
    def get_line_form(self):
        for index in range(2*self.__Nmax + 1):
            freq = self.__line.wavelength + self.__df*self.__Nmax - self.__df*index
            self.__amp_curve_freq.append(freq)
            self.__amp_curve_ampl.append(self.__amp_coeff*lineform(freq,
                            self.__line.wavelength, self.__half_line_width))
        pass

if __name__ == "__main__":
    line10P20 = wl.Line(10,'P', 20, unit='hz' )
    gas_mixture1 = GasMixture(temp=400)
    gas_volume1 = GasVolume()
    branch_10P = wl.Branch(10, 'P', 2, 58, unit='hz')
    for line in branch_10P:
        amplification = AmplificationCurve(line, gas_mixture = gas_mixture1, gas_volume = gas_volume1)
        print(line.name, ' - ',amplification.amplif)


    
