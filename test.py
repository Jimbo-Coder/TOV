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



gamma = 2 
K = 1

n = 1/(gamma-1);  # == 1


# Reading and creating 4 Layer Polytrope block. 
#  
def fourlayer():
   reader = csv.reader(open("/data/data_maxwell/TOV/peos_parameter.dat_SLy"))
   polyN = 4
   rhob =[]; gammas = []; presbs = [];counter = 1;
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

   #rhob = np.flip(rhob); gammas = np.flip(gammas);

   rhob = rhob * dsc *(100**3)/(1000); 
   rhozero = rhozero * dsc*(100**3)/(1000); preszero = preszero * psc*(100/1000);

   gammas = gammas[0:polyN];

   posi = np.digitize(rhozero, rhob)

   kappas = np.zeros(polyN);
   kappas[posi-1] = (preszero)/(rhozero**(gammas[posi-1]))

   for i in range(polyN):
      try:
         kappas[posi +i] = kappas[posi + i - 1] * (rhob[posi+i]**(gammas[posi+i-1]-gammas[posi+i]))
      except:
         break

   for i in range(polyN):
      if posi-2-i >=0:
         try:
               kappas[posi -2 - i] = kappas[posi -1-i] * (rhob[posi-1-i]**(gammas[posi-1-i]-gammas[posi-2-i]))
         except:
               break
      else:
         break

      
   deltas = []
   for j in range(polyN-1):
      d = kappas[j]*(rhob[j+1]**[gammas[j]]) - kappas[j+1]*(rhob[j+1]**[gammas[j+1]])
      deltas = np.append(deltas,d)

   resP = kappas * rhob[0:polyN]**gammas;
   resP = np.append(resP, kappas[polyN-1]*rhob[polyN]**gammas[polyN-1])
   print(deltas/resP[1:polyN])

   return(kappas, gammas, rhob, resP)

kappas, gammas, rhob, resP = fourlayer()

def Pres(b):
   #Basic 1 Layer polytrope
   #g = K* (b**gamma)

   # 4 Layer polytrope digitize the rest mass density to the bins rhob
    ilayer = np.digitize(b,rhob)
    ilayer = np.clip(ilayer, 1,4)
    g = kappas[ilayer-1] * (b **(gammas[ilayer-1]))
    return(g)