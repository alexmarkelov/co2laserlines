# Prpogramm to calculate wavelength of emission lines
# of CO2 molecule[cm-1] in transition 10mkm, 9mkm


import molconst12c16o2
import numpy as np


CO210MKM = (960.9586, 0.77733033, -3.047548082e-3, -4.95971792e-7,-1.8082522e-8, 6.24e-13, -1.74e-13)
CO29MKM = (1.0637346e3, 0.7776237, -3.340875873e-3, -5.8038975e-7, 2.4126352e-8, 7.5e-13, -2.16e-13)



class CO2LaserLine():

    
    MKM10 = molconst12c16o2.CO210MKM
    MKM9 = molconst12c16o2.CO29MKM
    
    def __init__(self, jnum, branch, central_wavelength):
        self.jnum = jnum
        self.branch = branch
        self.cen_wav = central_wavelength
        self.wavelength = 0
        self.getwavelength()
    
    def getwavelength(self):
        if self.cen_wav == 10:
            if self.branch == "P":
                self.wavelength = self.wavelengthcalc(-self.jnum, CO2LaserLine.MKM10)
                pass
            elif self.branch == "R":
                self.wavelength = self.wavelengthcalc(self.jnum+1, CO2LaserLine.MKM10)
                pass
        elif():
            if self.branch == "P":
                self.wavelength = self.wavelengthcalc(-self.jnum, CO2LaserLine.MKM9)
                pass
            elif self.branch == "R":
                self.wavelength = self.wavelenghtcalc(self.jnum+1, CO2LaserLine.MKM9)
                pass
        
    
    def wavelengthcalc(self, m, constants):
        wavelength = 0
        for member in enumerate(constants):
            wavelength = wavelength + m**member[0]*member[1]
        return wavelength


line10P20 = CO2LaserLine(20, "P", 10)
print(line10P20.wavelength)
 
 
