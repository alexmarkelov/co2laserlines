import matplotlib.pyplot as plt
import amplif_lines as amp
import wavelength as wl

gas_mixture1 = amp.GasMixture(temp=400)
gas_volume1 = amp.GasVolume()
branch_10P = wl.Branch(10, 'P', 2, 58, unit='hz')
x_line = list()
y_line = list()
fig, ax = plt.subplots()
for line in branch_10P:
    amplification = amp.AmplificationCurve(line, gas_mixture=gas_mixture1,
                                           gas_volume=gas_volume1)
    ax.plot(amplification.amp_curve_freq, amplification.amp_curve_ampl)
plt.show()
