import numpy as np
import scipy as sci
import matplotlib.pyplot as plt


G = 1; c = 1; 

gamma = 2 
K = 1

n = 1/(gamma-1)

def Pres(b):
   g = K* (b**gamma)
   return(g)

def density(z):
   g = (z/K)**(1/gamma)
   return(g)

def gradm(rhot, r):
   g = 4*np.pi*(r**2)*rhot
   return(g)

def edensity(rho, P):
   h = rho + P / (gamma-1)
   return(h)


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

dr = 0.0001; rc = 0.0001; rmax = 5;


def createstar(x, meth):

   rhoc = x; Pc = Pres(rhoc);  
   rhoct = edensity(rhoc,Pc)

   mc = rhoct*4*np.pi*(rc**2)
   mc0 = rhoc*4*np.pi*(rc**2)

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


rho0test = np.arange(0, .8, 0.01); Pctest = Pres(rho0test);

rhottest = edensity(rho0test, Pctest);

resultmass = []; resultmass0 = [];
for j in rho0test:
   mass, mass0,mfi,Pfi,gf = createstar(j, "Euler"); 
   resultmass = np.append(resultmass, mass)
   resultmass0 = np.append(resultmass0,mass0)


rhotvis = 1.47E15; rhotvis = rhotvis * (6.67E-11)*((3e8)**(-2)) * (100**(3));
rhotvis =0.42;
masst,masst0, mft, Pft,k = createstar(rhotvis, "Euler")


rn = (K**(-n/2))*rf;
rho0n = (K**(n))*rho0test;
rhotn = (K**(n))*rhottest;
mn = (K**(-n/2))*resultmass; masstn = (K**(-n/2))*masst; mftn = (K**(-n/2))*mft;
mn0 = (K**(-n/2))*resultmass0; masstn0 = (K**(-n/2))*masst;
rnf = rn[k];

ratio = masst/rnf;


peaki = np.argmax(mn)

print(f"Turning point found at rhoct = {rhotn[peaki]}, M = {mn[peaki]}, M_0 = {mn0[peaki]}")
print(f"For rho_c={rhotvis}, compactness is {ratio} ")

plt.figure()
plt.title("TOV M_Stars vs rho_c")
plt.ylabel("M_Star")
plt.xlabel("rhoc")
plt.axvline(x=rhotvis, c='k')
plt.plot(rhotn, mn)   
plt.plot(rhotn,mn0, linestyle="--")
plt.savefig("TOV M_Stars vs rho_c.pdf",dpi=300)


#sols = sol.sol(rs)
#rhos = sols[0]; ms = sols[1]



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
plt.plot(rn[0:k], mftn[0:k])
plt.savefig("TOV m vs r.pdf",dpi=300)

plt.figure()
plt.title("TOV P vs r")
plt.ylabel("Pressure")
plt.xlabel("Radius") 
plt.plot(rn[0:k], Pft[0:k])
#plt.show()
plt.savefig("TOV P vs r.pdf",dpi=300)
