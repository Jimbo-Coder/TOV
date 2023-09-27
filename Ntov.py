import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import csv


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

gamma = 2 
K = 1

n = 1/(gamma-1);  # == 1


# Reading and creating 4 Layer Polytrope block. 
#  
def fourlayer():
   reader = csv.reader(open("peos_parameter.dat_SLy"))
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
   rhozero = rhozero * dsc*(100**3)/(1000); preszero = preszero * psc;

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

      
   for j in range(polyN-1):
      d = kappas[j]*(rhob[j+1]**[gammas[j]]) - kappas[j+1]*(rhob[j+1]**[gammas[j+1]])

   resP = kappas * rhob[0:polyN]**gammas;
   print(d/resP[1:polyN])

   return(kappas, gammas, rhob)

kappas, gammas, rhob = fourlayer()
# End Read Polytrope Block



# Currently rhob = [1.62820830e+00, 1.62820830e-03, 8.16036834e-04, 2.38076618e-04, 1.62820830e-18]
# 1.62820830 limit
def Pres(b):
   #Basic 1 Layer polytrope
   #g = K* (b**gamma)

   # 4 Layer polytrope digitize the rest mass density to the bins rhob

   ilayer = np.digitize(b,rhob)

   ilayer = np.clip(ilayer, 1,4)
   g = kappas[ilayer-1] * (b **(gammas[ilayer-1]))
   
   return(g)





def density(z):
   g = (z/K)**(1/gamma)
   return(g)

def gradm(rhot, r):
   g = 4*np.pi*(r**2)*rhot
   return(g)

def gradm0(rho, m, r):
   h = 4*np.pi*(r**2)*rho*((1- (2*m)/r)**(-1/2))
   return(h)

def gradp(rho,rhot, m, r):
   h = ((-rhot*m)/(r**2)) * (1 + Pres(rho)/rhot) * (1 + (4 * np.pi * Pres(rho)*(r**3))/m)
   f = h * ((( 1 - ((2* m) / r)))**(-1))
   return(f)



def schwarz(M, R):
   g = 0.5* np.log(1 - ((2*M)/R))
   return(g)

def edensity(rho, P):
   h = rho + P / (gamma-1)
   return(h)



dr = 0.0001; rc = 0.000000001; rmax = 100;

def createstar(x, meth):

   rhoc = x; Pc = Pres(rhoc);  
   rhoct = edensity(rhoc,Pc)

   mc = (4/3)*np.pi*(rc**3)*rhoc
   mc0 = (4/3)*np.pi*(rc**3)*rhoc

   rf = np.arange(rc, rmax, dr, dtype = np.float64)

   mf = np.zeros(len(rf));Pf = np.zeros(len(rf));rhof = np.zeros(len(rf));
   rhoft = np.zeros(len(rf)); mf0 = np.zeros(len(rf));

   mf[0]+=mc; Pf[0] += Pc; rhof[0] += rhoc; rhoft[0] += rhoct; mf0[0] += mc0



   for i in range(len(rf)-1):

      if (rhof[i] ==0) or (mf[i] == 0)  or (rf[i]== 0) or (rhoft[i]==0) or (mf0[i]==0):
         mres = 0; phif0 = 0; 
         break
      if meth == "Euler":

         dm = gradm(rhoft[i],rf[i])
         dm0 = gradm0(rhof[i],mf[i],rf[i])
         dP = gradp(rhof[i],rhoft[i],mf[i],rf[i])

         Pf[i+1] = Pf[i] + dP * dr
         mf[i+1] = mf[i] + dm * dr
         mf0[i+1] = mf0[i] + dm0 * dr

         if Pf[i+1] <= 0:
            break
      
         rhof[i+1] = density(Pf[i+1])
         rhoft[i+1] = edensity(rhof[i+1],Pf[i+1])

      elif meth == "RK4":
         
         dm = gradm(rhoft[i],rf[i])
         dm0 = gradm0(rhof[i],mf[i],rf[i])
         dP = gradp(rhof[i],rhoft[i],mf[i],rf[i])

         Pf[i+1] = Pf[i] + dP * dr
         mf[i+1] = mf[i] + dm * dr
         mf0[i+1] = mf0[i] + dm0 * dr

         if Pf[i+1] <= 0:
            break

         rhof[i+1] = density(Pf[i+1])
         rhoft[i+1] = edensity(rhof[i+1],Pf[i+1])

   mres = mf[i]
   mres0 = mf0[i]
   return(mres,mres0, mf,Pf,i)



rf = np.arange(rc, rmax, dr, dtype = np.float64)

#1.791287
rho0test = np.arange(0, 1.79, 0.01); Pctest = Pres(rho0test);

rhottest = edensity(rho0test, Pctest);

resultmass = []; resultmass0 = [];
for j in rho0test:
   mass, mass0,mfi,Pfi,gf = createstar(j, "Euler"); 
   resultmass = np.append(resultmass, mass)
   resultmass0 = np.append(resultmass0,mass0)


rhotvis = 1.47E15; rhotvis = rhotvis*dsc *(1/1000)*(100**3);# 0.0023934662036213987
masst,masst0, mft, Pft,k = createstar(rhotvis, "Euler")




rnf = rf[k];

ratio = masst/rnf;


peaki = np.argmax(resultmass)

print(f"Turning point found at rhoct = {rhottest[peaki]}, M = {resultmass[peaki]}, M_0 = {resultmass0[peaki]}")
print(f"For rho_c={rhotvis}, Radius is {rf[k]}, M_t  = {masst}, M_0 = {masst0},compactness is {ratio} ")

counter = 1;
reader = csv.reader(open("tov.dat_orig"))
pa = []; ma = []; mb = []; mc = []; p0 = []
for row in reader:
   if counter <=7:
      counter +=1
      continue
   else:
      b = row[0].split()
      a = b[0], b[2],b[3],b[4],b[6]
      pa = np.append(pa, a[0])
      ma = np.append(ma, a[1])
      mb = np.append(mb, a[2])
      mc = np.append(mc, a[3])
      p0 = np.append(p0, a[4])

   counter +=1

ma = np.array(ma, dtype = np.float64);
mb = np.array(mb, dtype = np.float64);
mc = np.array(mc, dtype = np.float64);
pa = np.array(pa, dtype = np.float64);
p0 = np.array(p0, dtype = np.float64);


plt.figure()
plt.title("TOV M_Stars vs rho_c")
plt.ylabel("M_Star")
plt.xlabel("rhoc")
plt.axvline(x=rhotvis, c='b')
plt.plot(rhottest, resultmass, label = "mM_t")   
plt.plot(rhottest,resultmass0, linestyle="--", label = "mM_0")

plt.scatter(pa, ma,label = "M_ADM", s = 0.65,c='k'); 
plt.scatter(pa, mb, label = "M_prop", s = 0.65,c='m'); 
plt.scatter(pa, mc, label = "M_0", s = 0.65,c='y');
plt.xlim(0.01,1.6)
plt.legend()
plt.savefig("TOV M_Stars vs rho_c.pdf",dpi=300)


# RESIDUAL PLOTTING BLOCK, Uncomment if want comparisons

# resmass = []; resmass0 = [];
# for j in p0:
#    mass, mass0,mfi,Pfi,gf = createstar(j, "Euler"); 
#    resmass = np.append(resmass, mass)
#    resmass0 = np.append(resmass0,mass0)



# plt.figure()
# plt.title("ADM/Total Mass Residual")
# plt.plot(pa, ((resmass - ma)/ma))
# plt.ylabel("deltaM")
# plt.xlabel("Central Total Energy Density")
# plt.savefig("ADM_MT resid.pdf",dpi=300)

# plt.figure()
# plt.title("Rest Mass resid")
# plt.plot(pa, ((resmass0 - mc)/mc))
# plt.ylabel("deltaM")
# plt.xlabel("Central Total Energy Density")
# plt.savefig("Rest Mass Resid.pdf",dpi=300)


##END RESIDUAL BLOCK





# plt.figure()
# plt.title("TOV rho vs r")
# plt.ylabel("Density")
# plt.xlabel("Radius")
# plt.plot(rf[0:i], rhof[0:i])
# plt.show()


plt.figure()
plt.title("TOV m vs r")
plt.ylabel("Mass")
plt.xlabel("Radius")
plt.plot(rf[0:k], mft[0:k])
plt.savefig("TOV m vs r.pdf",dpi=300)

plt.figure()
plt.title("TOV P vs r")
plt.ylabel("Pressure")
plt.xlabel("Radius") 
plt.plot(rf[0:k], Pft[0:k])
plt.savefig("TOV P vs r.pdf",dpi=300)
