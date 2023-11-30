import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import pandas as pd

G = 6.67430e-11; # 6.67E-11 m^3 / kg*s^2 = 1
c = 2.99792458e8; # 3e8 m = 1 second
Msun = 1.989e30; #2E30 Kg = 1, call this Msun
lsc = G*Msun / (c**2) # so GMsun/c^2 = 1,482.222 m = 1, call this lsc
tsc = G*Msun / (c**3) # then GMsun/c^3 = 4.94e-6 s = 1, call this tsc
#energy is M*L^2/T^2 , mult by tsc^2 / (lsc^2 Msun) and its pure number
psc = lsc*(tsc**2)/Msun #Pressure = energydensity, multi by  lsc*tsc^2 / (Msun) and its pure, factor psc
dsc = (lsc**3)/Msun #For normal mass density
#Scalings are Msun, lsc, tsc, psc
# G = C = M_Sun  = 1 UNITS
function readtabeos():

    path = "/data/data_maxwell/TOV/teos_parameter.dat_Togashi"
    tab = pd.read_csv(path,delim_whitespace='True', header = None,skiprows = 1)
    tabef = tab[0]; tabpf = tab[1]
    return(tabef, tabpf)

tabe, tabp = readtabeos()

np.digitize()