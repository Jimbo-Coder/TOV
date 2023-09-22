import numpy as np
import os
import csv
import matplotlib.pyplot as plt


G = 6.67e-11; # 6.67E-11 m^3 / kg*s^2 = 1
c = 3e8; # 3e8 m = 1 second
Msun = 2e30; #2E30 Kg = 1, call this Msun
lsc = G*Msun / (c**2) # so GMsun/c^2 = 1,482.222 m = 1, call this lsc
tsc = G*Msun / (c**3) # then GMsun/c^3 = 4.94e-6 s = 1, call this tsc
#energy is M*L^2/T^2 , mult by tsc^2 / (lsc^2 Msun) and its pure number
psc = lsc*(tsc**2)/Msun #Pressure = energydensity, multi by  lsc*tsc^2 / (Msun) and its pure, factor psc
dsc = (lsc**3)/Msun #For normal mass density
#Scalings are Msun, lsc, tsc, psc
# G = C = M_Sun  = 1 UNITS

reader = csv.reader(open("peos_parameter.dat_SLy"))
polyN = 4
rhob =[]; gammas = []; presbs = [];counter = 1
for row in reader:
    if counter == 8:
        break
    a = row[0]
    a = a.split()
    if counter == 1:
        counter +=1
    elif counter ==2:
        rhozero = np.array(a[0], dtype = np.float64)
        preszero = np.array(a[1], dtype = np.float64)
        counter+=1
    elif counter>=2:
        rhob = np.append(rhob, a[0])
        gammas = np.append(gammas,a[1])
        counter+=1
    

rhob = np.array(rhob, dtype = np.float64);
gammas = np.array(gammas, dtype = np.float64);

rhob = np.flip(rhob); gammas = np.flip(gammas);

rhob = rhob * dsc; rhozero = rhozero * dsc; preszero = preszero * psc;

rhob = rhob[1:5]; gammas = gammas[1:5]

posi = np.digitize(rhozero, rhob)

kappas = np.zeros(4);

for i in range(int(polyN/2)):
    if i==0:
        kappas[posi + i] = (preszero)/rhozero*(gammas[posi + i])
        kappas[posi -1-i] = (preszero)/rhozero*(gammas[posi -1 -i])
    else:
        kappas[posi + i] = kappas[posi + i - 1] * (rhob[posi + i - 1]**(gammas[posi + i - 1]-gammas[posi + i]))
        kappas[posi -1-i] = kappas[posi -i] * (rhob[posi -i]**(gammas[posi -i]-gammas[posi - i-1]))

    
for j in range(3):
    d = kappas[j]*(rhob[j+1]**[gammas[j]]) - kappas[j+1]*(rhob[j+1]**[gammas[j+1]])
    print(d)

g = kappas * (rhob**gammas)
plt.figure()
plt.plot(g)