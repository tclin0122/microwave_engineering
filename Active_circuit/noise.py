import numpy as np
import cmath
import scipy.constants as sy

B = 20 * np.power(10,6)
k = sy.Boltzmann
T = 290
Pn = k*T*B
Pn_dB = 10*np.log10(Pn*1000)
print("Available noise power = ",Pn_dB,"dBm")

def cascadeLNA(Pn,Ncascade,G_dB,G_prev,NF,F_tot,stage):
    Gav = np.power(10,(G_dB/10))
    F = np.power(10,(NF/10))
    Namp = Gav*(F-1)*Pn
    N = Ncascade*Gav + Namp
    N_dB = 10 * np.log10(N*1000)
    if(stage == 1): 
        F_tot = F
    else:
        F_tot = F_tot + (F - 1)/G_prev
    NF_tot = 10 * np.log10(F_tot)
    return N, Gav, N_dB, F_tot, NF_tot

Gav_dB_1 = 12
NF_1 = 2
Gav_dB_2 = 30
NF_2 = 7

(N1, Gav_1, N1_dB, F_tot,  NF_tot) = cascadeLNA(Pn, Pn, Gav_dB_1, 0, NF_1,0,1)
(N2, Gav_2, N2_dB, F_tot,  NF_tot) = cascadeLNA(Pn, N1, Gav_dB_2, Gav_1,NF_2,F_tot,2)

print("Available noise power at the input of the buffer amplifier = ", N1_dB, "dBm")
print("Total noise figure = ", NF_tot,"dB")

Sin = -70 #dBm
SNR = Sin - Pn_dB - NF_tot
print("Signal-to-noise ratio at the output of the buffer amplifier.", SNR, "dB")


print("")
F = np.power(10,(0.5/10))
F_min = np.power(10,(0.15/10))
Z_opt = complex(72,27)
gamma_opt = (Z_opt - 50)/(Z_opt + 50)
Rn = 2

N = (F - F_min)/(4*Rn/50)*np.power(abs(1+gamma_opt),2)
O_N = (gamma_opt)/(1+N)
R_N = 1/(1+N)*np.sqrt(N*N+N*(1-np.power(abs(gamma_opt),2)))


print("Center (absolute value)  for a noise figure of 0.5 dB = ", abs(O_N))
print("Center (phase in degrees) for a noise figure of 0.5 dB = ", np.rad2deg(cmath.phase(O_N)) )
print("Radius of the noise circle = ",R_N)

'''
F = np.power(10,(1/10))
F_min = np.power(10,(0.57/10))
Z_opt = complex(100,5.2)
gamma_opt = (Z_opt - 50)/(Z_opt + 50)
print(gamma_opt)
Rn = 6

N = (F - F_min)/(4*Rn/50)*np.power(abs(1+gamma_opt),2)
print(N)
O_N = (gamma_opt)/(1+N)
R_N = 1/(1+N)*np.sqrt(N*N+N*(1-np.power(abs(gamma_opt),2)))

print(O_N)
print(R_N)
'''