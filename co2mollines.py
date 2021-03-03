# Prpogramm to calculate wavelength of emission lines
# of CO2 molecule[cm-1] in transition 10mkm, 9mkm


import molconst12c16o2
import numpy as np


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
        elif self.cen_wav == 9:
            if self.branch == "P":
                self.wavelength = self.wavelengthcalc(-self.jnum, CO2LaserLine.MKM9)
                pass
            elif self.branch == "R":
                self.wavelength = self.wavelengthcalc(self.jnum+1, CO2LaserLine.MKM9)
                pass
        
    
    def wavelengthcalc(self, m, constants):
        wavelength = 0
        for member in enumerate(constants):
            wavelength = wavelength + m**member[0]*member[1]
        return wavelength


def printalllines():
    branches = ("P", "R")
    centers = (10, 9)
    for branch in branches:
        for center in centers:
            for line in range(1,25):
                A = CO2LaserLine(2*line, branch, center)
                print(A.cen_wav, A.branch, A.jnum,' wavelength[cm-1] - ', A.wavelength)
        

if __name__ == '__main__':
    printalllines()


 
 
