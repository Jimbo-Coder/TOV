import numpy as np
import scipy as sci
import matplotlib.pyplot as plt


G = 1; c = 1; 

γ = 2
K = 1

n = 1/(γ-1)

def Pres(ρt):
   g = K* (ρt**γ)
   return(g)

def density(z):
   g = (z/K)**(1/γ)
   return(g)

def gradm(ρt, r):
   g = 4*np.pi*(r**2)*ρt
   return(g)

def edensity(ρ, P):
   h = ρ + P / (γ-1)
   return(h)



def gradp(ρ,ρt, m, r):
   h = ((-ρt*m)/(r**2)) * (1 + Pres(ρ)/ρt) * (1 + (4 * np.pi * Pres(ρ)*(r**3))/m)
   f = h * ((( 1 - ((2* m) / r)))**(-1))
   return(f)

def schwarz(M, R):
   g = 0.5* np.log(1 - ((2*M)/R))
   return(g)

dr = 0.0001; rc = 0.0001; rmax = 5;


def createstar(x):

   ρc = x; Pc = Pres(ρc);  
   ρct = edensity(ρc,Pc)

   mc = ρct*4*np.pi*(rc**2)
   mc0 = ρc*4*np.pi*(rc**2)

   rf = np.arange(rc, rmax, dr, dtype = np.float64)

   mf = np.zeros(len(rf));Pf = np.zeros(len(rf));ρf = np.zeros(len(rf));
   ρft = np.zeros(len(rf)); mf0 = np.zeros(len(rf));

   mf[0]+=mc; Pf[0] += Pc; ρf[0] += ρc; ρft[0] += ρct; mf0[0] += mc0



   for i in range(len(rf)-1):

      if (ρf[i] ==0) or (mf[i] == 0)  or (rf[i]== 0) or (ρft[i]==0) or (mf0[i]==0):
         mres = 0; phif0 = 0; mres0 = 0;
         break

      dm = gradm(ρft[i],rf[i])
      dm0 = gradm(ρf[i],rf[i])
      dP = gradp(ρf[i],ρft[i],mf[i],rf[i])

      Pf[i+1] = Pf[i] + dP * dr
      mf[i+1] = mf[i] + dm * dr
      mf0[i+1] = mf0[i] + dm0 * dr

      if Pf[i+1] <= 0:
         break
      ρf[i+1] = density(Pf[i+1])
      ρft[i+1] = edensity(ρf[i+1],Pf[i+1])

   mres = mf[i]
   mres0 = mf0[i]
   return(mres,mres0, mf, Pf,i)



rf = np.arange(rc, rmax, dr, dtype = np.float64)


ρ0test = np.arange(0, 1.3, 0.01); Pctest = Pres(ρ0test);

ρttest = edensity(ρ0test, Pctest);

resultmass = []; resultmass0 = [];
for j in ρ0test:
   mass, mass0 = createstar(j)[0:2]; 
   resultmass = np.append(resultmass, mass)
   resultmass0 = np.append(resultmass0,mass0)


ρtvis = 1.47E15; ρtvis = ρtvis * (6.67E-11)*((3e8)**(-2)) * (100**(3));
ρtvis =0.42;
masst,masst0, mft, Pft,k = createstar(ρtvis)


rn = (K**(-n/2))*rf;
ρ0n = (K**(n))*ρ0test;
ρtn = (K**(n))*ρttest;
mn = (K**(-n/2))*resultmass; masstn = (K**(-n/2))*masst; mftn = (K**(-n/2))*mft;
mn0 = (K**(-n/2))*resultmass0; masstn0 = (K**(-n/2))*masst;
rnf = rn[k];
ratio = masst/rnf;


peaki = np.argmax(mn)

print(f"Turning point found at ρct = {ρtn[peaki]}, M = {mn[peaki]}, M_0 = {mn0[peaki]}")
print(f"For ρ={ρtvis}, compactness is {ratio} ")

plt.figure()
plt.title("TOV M_Stars vs ρ_c")
plt.ylabel("M_Star")
plt.xlabel("ρc")
plt.axvline(x=ρtvis, c='k')
plt.plot(ρtn, mn)   
plt.plot(ρtn,mn0, linestyle="--")
plt.show()


#sols = sol.sol(rs)
#ρs = sols[0]; ms = sols[1]



# plt.figure()
# plt.title("TOV ρ vs r")
# plt.ylabel("Density")
# plt.xlabel("Radius")
# plt.plot(rf[0:i], ρf[0:i])
# plt.show()

plt.figure()
plt.title("TOV m vs r")
plt.ylabel("Mass")
plt.xlabel("Radius")
plt.plot(rn[0:k], mftn[0:k])
plt.show()

plt.figure()
plt.title("TOV P vs r")
plt.ylabel("Pressure")
plt.xlabel("Radius")
plt.plot(rn[0:k], Pft[0:k])
plt.show()
