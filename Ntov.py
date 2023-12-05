import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import csv
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
# End Read Polytrope Block

#Read Tabulated Equation of State block
def readtabeos():
    pathtab = "/data/data_maxwell/TOV/teos_parameter.dat_Togashi"
    tab = pd.read_csv(pathtab,delim_whitespace='True', header = None,skiprows = 1)
    tabef = tab[0]; tabpf = tab[1]; tabp0f = tab[3];
    return(tabef,tabp0f, tabpf)
#Call
tabe, tabp0, tabp = readtabeos()



# Currently 4layer rhob = [1.62820830e+00, 1.62820830e-03, 8.16036834e-04, 2.38076618e-04, 1.62820830e-18]
# 1.62820830 limit
#Def Pressure function. 4layer polytrope or tabulated equation of state.
#Pressure as a function of rest-mass density, and the equation of state
def Pres(b, eosm):
   #Basic 1 Layer polytrope
   #g = K* (b**gamma)

   # 4 Layer polytrope digitize the rest mass density to the bins rhob
   if eosm == "4layer":
      ilayer = np.digitize(b,rhob)
      ilayer = np.clip(ilayer, 1,4)
      g = kappas[ilayer-1] * (b **(gammas[ilayer-1]))
      return(g)

   elif eosm == "tab":
      ilayer = np.digitize(b,tabp0)
      ilayer = np.clip(ilayer, 1,5499)
      tp = tabp[ilayer-1]; tp2 = tabp[ilayer];
      tpi = np.interp(b, [ilayer,ilayer-1],[tp,tp2])
      return(tpi)



#Affix to Onelayer
#kappas = np.array([1,1,1,1]); gammas = np.array([2,2,2,2])

#Function that returns the Rest-mass density, given the Pressure and eos
def density(z, eosm):
   if eosm == "4layer":
      ilayer = np.digitize(z,resP)
      ilayer = np.clip(ilayer, 1,4)
      g = (z/kappas[ilayer-1])**(1/gammas[ilayer-1])
      return(g)

   elif eosm == "tab":
      ilayer = np.digitize(z,tabp)
      ilayer = np.clip(ilayer, 1,5499)
      tp = tabp0[ilayer-1]; tp2 = tabp0[ilayer];
      tpi = np.interp(z, [ilayer,ilayer-1],[tp,tp2])
      return(tpi)

#Function that returns the total energy density, given p0, P, and eos
def edensity(rho, P,eosm):
   if eosm == "4layer":
      ilayer = np.digitize(rho,rhob)
      ilayer = np.clip(ilayer, 1,4)
      h = rho + P / (gammas[ilayer-1]-1)
      return(h)

   elif eosm == "tab":
      ilayer = np.digitize(rho,tabp0)
      ilayer = np.clip(ilayer, 1,5499)
      tp = tabe[ilayer-1]; tp2 = tabe[ilayer];
      tpi = np.interp(rho, [ilayer,ilayer-1],[tp,tp2])
      return(tpi)
   
#dm/dr = 4pi*r^2*rho
def gradm(rhot, r):
   g = 4*np.pi*(r**2)*rhot
   return(g)

#dm0/dr = 4pi*r^2*rho*g_00^0.5
def gradm0(rho, m, r):
   h = 4*np.pi*(r**2)*rho*((1- (2*m)/r)**(-1/2))
   return(h)

#TOV dP/dr equation, needs both densities, mass, radius, and eos
def gradp(rho,rhot, m, r, a):
   h = ((-rhot*m)/(r**2)) * (1 + Pres(rho,a)/rhot) * (1 + (4 * np.pi * Pres(rho,a)*(r**3))/m)
   f = h * ((( 1 - ((2* m) / r)))**(-1))
   return(f)

#Not used
def schwarz(M, R):
   g = 0.5* np.log(1 - ((2*M)/R))
   return(g)



def createstar(x, meth,metheos):

   rhoc = x; Pc = Pres(rhoc, metheos);  
   rhoct = edensity(rhoc,Pc,metheos)

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
         dP = gradp(rhof[i],rhoft[i],mf[i],rf[i],metheos)

         Pf[i+1] = Pf[i] + dP * dr
         mf[i+1] = mf[i] + dm * dr
         mf0[i+1] = mf0[i] + dm0 * dr

         if Pf[i+1] <= 0:  #If pressure is negative, break
            break
      
         rhof[i+1] = density(Pf[i+1],metheos)
         rhoft[i+1] = edensity(rhof[i+1],Pf[i+1],metheos)



      elif meth == "RK4":
         dm = np.zeros(4);dm0 = np.zeros(4); dP = np.zeros(4);

         dm[0] = gradm(rhoft[i],rf[i])*dr
         dm0[0] = gradm0(rhof[i],mf[i],rf[i])*dr
         dP[0] = gradp(rhof[i],rhoft[i],mf[i],rf[i],metheos)*dr

         if Pf[i] + dP[0]/2 <= 0:
            break

         rho1t = density(Pf[i] + dP[0]/2,metheos)
         rho1tt = edensity(rho1t,Pf[i]+ dP[0]/2,metheos)

         dP[1] = gradp(rho1t,rho1tt,mf[i] + dm[0]/2,rf[i] + dr/2,metheos)*dr
         dm[1] = gradm(rho1tt,rf[i] + dr/2)*dr
         dm0[1] = gradm0(rho1t,mf[i]+ dm[0]/2,rf[i] + dr/2)*dr
         
         if Pf[i] + dP[1]/2 <= 0:
            break

         rho2t = density(Pf[i] +  dP[1]/2,metheos)
         rho2tt = edensity(rho1t,Pf[i]+ dP[1]/2,metheos)

         dP[2] = gradp(rho2t,rho2tt,mf[i] + dm[1]/2,rf[i] + dr/2,metheos)*dr   
         dm[2] = gradm(rho2tt,rf[i] + dr/2)*dr
         dm0[2] = gradm0(rho2t,mf[i]+ dm[1]/2,rf[i] + dr/2)*dr
         
         if Pf[i] + dP[2] <= 0:
            break

         rho3t = density(Pf[i] +  dP[2],metheos)
         rho3tt = edensity(rho1t, Pf[i]+ dP[2],metheos)

         dP[3] = gradp(rho3t,rho3tt,mf[i] + dm[2],rf[i] + dr,metheos)*dr   
         dm[3] = gradm(rho3tt,rf[i] + dr)*dr
         dm0[3] = gradm0(rho3t,mf[i]+ dm[2],rf[i] + dr)*dr

         
         Pf[i+1] = Pf[i] + (dP[0] + 2*dP[1] + 2*dP[2] + dP[3])/6
         mf[i+1] = mf[i] + (dm[0] + 2*dm[1] + 2*dm[2] + dm[3])/6
         mf0[i+1] = mf0[i] + (dm0[0] + 2*dm0[1] + 2*dm0[2] + dm0[3])/6 

         if Pf[i+1] <= 0:
            break

         rhof[i+1] = density(Pf[i+1],metheos)
         rhoft[i+1] = edensity(rhof[i+1],Pf[i+1],metheos)






   mres = mf[i]
   mres0 = mf0[i]
   return(mres,mres0, mf,Pf,i)

dr = 0.01; rc = 0.0000001; rmax = 300;
rf = np.arange(rc, rmax, dr, dtype = np.float64)

#1.791287
#rho0test = np.arange(0, 1.78, 0.01); Pctest = Pres(rho0test);

#A few simple tests to ensure results are reasonable
rhotvis = 1.47E15; rhotvis = rhotvis*dsc *(1/1000)*(100**3);# 0.0023934662036213987
rhotvis = 1.42E15; rhotvis = rhotvis*dsc *(1/1000)*(100**3); # 0.00230044

#Run a single simulation on the test values
masst,masst0, mft, Pft,k = createstar(rhotvis, "RK4","4layer")

rnf = rf[k];

ratio = masst/rnf;

#Read 1layer polytrope data
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
#End 1layer polytrope data


#Read 4 layer polytrope data
path = "ovphy_plot_SLy_EOS_00.dat"

b = pd.read_csv(path,header=None,delim_whitespace=True)
rhoant = b[2]; m0ant = b[8]; madant = b[10];rsant = b[13];
#End 4 layer polytrope data


#Iterateover the given densities, and create a star for each one
rho0test = np.arange(0, rhoant[len(rhoant)-1], 0.001); 
rho0test = rhoant
Pctest = Pres(rho0test, "4layer"); rhottest = edensity(rho0test, Pctest,"4layer");

resultmass = []; resultmass0 = []; resultmassrk=[]; resultmassrk0=[];
resrad = []; resradrk = [];resultmassrt = []; resultmassrt0 = []; resradrt=[];
resultmasset = []; resultmasset0 = []; resradet=[];
#Run the sim, 4 options of eos and Integration method
for j in rho0test:
   mass, mass0,mfi,Pfi,gf = createstar(j, "Euler","4layer"); 
   massrk,massrk0,mfi1,Pfi1,gf1 = createstar(j, "RK4","4layer");
   massrt,massrt0,mfi2,Pfi2,gft2 = createstar(j, "RK4","tab");
   masset, masset0, mfi3, Pfi3, gft3 = createstar(j, "Euler","tab");

   resultmass = np.append(resultmass, mass)
   resultmass0 = np.append(resultmass0,mass0)
   
   resultmassrt = np.append(resultmassrt, massrt)
   resultmassrt0 = np.append(resultmassrt0,massrt0)

   resultmasset = np.append(resultmasset, masset)
   resultmasset0 = np.append(resultmasset0,masset0)

   resultmassrk = np.append(resultmassrk, massrk) 
   resultmassrk0 = np.append(resultmassrk0,massrk0)

   resradrt = np.append(resradrt, rf[gft2])
   resrad = np.append(resrad, rf[gf])
   resradrk = np.append(resradrk, rf[gf1])
   resradet = np.append(resradet, rf[gft3])

   print(f'{j} star done')


peaki = np.argmax(resultmass)
peaki1 = np.argmax(resultmassrk)
#Find/print turning point
print(f"Turning point found at rhoct = {rhoant[peaki]}, M = {resultmass[peaki]}, M_0 = {resultmass0[peaki]}")
print(f"RK4 Turning point found at rhoct = {rhoant[peaki1]}, M = {resultmassrk[peaki1]}, M_0 = {resultmassrk0[peaki1]}")
print(f"For rho_c={rhotvis}, Radius is {rf[k]}, M_t  = {masst}, M_0 = {masst0}, M/R is {ratio} ")

#plot Rest mass residual
plt.figure()
plt.title("4Layer M_0 Resid")
plt.xlabel("rho_c")
plt.ylabel("rel err.")
plt.plot(rhoant, (resultmass0 - m0ant)/m0ant, label = "Euler")
plt.plot(rhoant, (resultmassrk0 - m0ant)/m0ant,label = "RK4")
plt.plot(rhoant, (resultmassrt0 - m0ant)/m0ant,label = "RK4tab")
plt.plot(rhoant, (resultmasset0 - m0ant)/m0ant,label = "Eultab")
plt.legend()
plt.savefig("4layer restm Resid")

#Plot ADM Mass residual
plt.figure()
plt.title("4Layer M_ADM Resid")
plt.xlabel("rho_c")
plt.ylabel("rel err.")
plt.plot(rhoant, (resultmass - madant)/madant, label = "Euler")
plt.plot(rhoant, (resultmassrk - madant)/madant, label = "RK4")
plt.plot(rhoant, (resultmassrt - madant)/madant, label = "RK4tab")
plt.plot(rhoant, (resultmasset - madant)/madant, label = "Eulertab")
plt.legend()
plt.savefig("4layer ADM Resid")

#Plot Radius residual
plt.figure()
plt.title("4Layer Radius Resid")
plt.xlabel("rho_c")
plt.ylabel("rel err.")
plt.plot(rhoant, (resrad - rsant)/rsant, label = "Euler")
plt.plot(rhoant, (resradrk - rsant)/rsant, label = "RK4")
plt.plot(rhoant, (resradrt - rsant)/rsant, label = "RK4tab")
plt.plot(rhoant, (resradet - rsant)/rsant, label = "Eulertab")
plt.legend()
plt.savefig("4layer radius Resid")

#Plot rk4/euler differences
plt.figure()
plt.title("Euler/RK4 Resid")
plt.xlabel("rho_c")
plt.ylabel("rel err.")
plt.plot(rhoant, (resultmassrk - resultmass)/resultmassrk, label = "M_ADM resid")
plt.plot(rhoant, (resultmassrk0 - resultmass0)/resultmassrk0, label = "M_0 resid")
plt.plot(rhoant, (resultmassrt - resultmass)/resultmassrt, label = "M_ADM Rk4Tab resid")
plt.plot(rhoant, (resultmassrt0 - resultmass0)/resultmassrk0, label = "M_0 Rk4Tab resid")
plt.plot(rhoant, (resultmasset - resultmassrk)/resultmassrk, label = "M_ADM EulTab resid")
plt.plot(rhoant, (resultmasset0 - resultmassrk)/resultmassrk, label = "M_ADM EulTab resid")
plt.legend()
plt.savefig("Euler_RK4 Resid")


#Main Plot of final stellar mass vs central density
plt.figure()
plt.title("TOV M_Stars vs rho_c")
plt.ylabel("M_Star")
plt.xlabel("ρ_c")
#plt.axvline(x=rhotvis, c='b')

#plt.plot(rho0test, resultmass, label = "EulerM_ADM")   
#plt.plot(rho0test,resultmass0, linestyle="--", label = "EulerM_0")

plt.plot(rho0test, resultmassrk, label = "Rk4M_ADM",c='b')   
plt.plot(rho0test,resultmassrk0, linestyle="--", label = "Rk4M_0",c='b')

plt.plot(rho0test, resultmassrt, label = "Rk4tabM_ADM",c='g')
plt.plot(rho0test, resultmassrt0, linestyle="--", label = "Rk4tabM_0",c='g')

plt.plot(rho0test, resultmasset, label = "EultabM_ADM",c='m')
plt.plot(rho0test, resultmasset0, linestyle="--", label = "EultabM_0",c='m')

plt.scatter(rhoant,m0ant, label = "M_0", s = 0.85,c='r')
plt.scatter(rhoant,madant, label = "M_ADM",s = 0.85,c='c')

#Comparison to tov.dat file

# plt.scatter(pa, ma,label = "M_ADM", s = 0.65,c='k'); 
# plt.scatter(pa, mb, label = "M_prop", s = 0.65,c='m'); 
# plt.scatter(pa, mc, label = "M_0", s = 0.65,c='y');
# #plt.xlim(0.01,1.6)
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
