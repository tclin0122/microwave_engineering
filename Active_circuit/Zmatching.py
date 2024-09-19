import numpy as np
import cmath
import scipy.constants as sy

## S parameters
Mag_S11 = 0.47
Ang_S11 = 173
Mag_S12 = 0
Ang_S12 = 0
Mag_S21 = 6.01
Ang_S21 = 55
Mag_S22 = 0.27
Ang_S22 = -75


## Convert S parameters to complex numbers
S11 = complex(Mag_S11 * np.cos(np.deg2rad(Ang_S11)),Mag_S11 * np.sin(np.deg2rad(Ang_S11)))
S21 = complex(Mag_S21 * np.cos(np.deg2rad(Ang_S21)),Mag_S21 * np.sin(np.deg2rad(Ang_S21)))
S12 = complex(Mag_S12 * np.cos(np.deg2rad(Ang_S12)),Mag_S12 * np.sin(np.deg2rad(Ang_S12)))
S22 = complex(Mag_S22 * np.cos(np.deg2rad(Ang_S22)),Mag_S22 * np.sin(np.deg2rad(Ang_S22)))

Z0 =50 #impedance
Zload = 50.00000000001
Zsource = 50.00000000001
ZL = (Zload-Z0) / (Zload+Z0)
ZS = (Zsource-Z0) / (Zsource+Z0)
refl_input = S11 + (S12*S21)/(1/ZL - S22)
#print("The absolute value of the input reflection coefficient = ",abs(refl_input))
refl_output = S22 + (S12*S21)/(1/ZS - S11)
#print("The absolute value of the output reflection coefficient = ",abs(refl_output))

## Gain calculation
GT = ((1 - np.power(abs(ZL),2))/(np.power(abs(1 - S22*ZL),2))) * np.power(abs(S21),2) * ((1 - np.power(abs(ZS),2))/(np.power(abs(1 - refl_input*ZS),2)))
GT_dB = 10*np.log10(GT)
Gav = ((1 - np.power(abs(ZS),2))/(np.power(abs(1 - S11*ZS),2))) * np.power(abs(S21),2) * (1/(1 - np.power(abs(refl_output),2)))
Gav_dB = 10*np.log10(Gav)
Gdel = ((1 - np.power(abs(ZL),2))/(np.power(abs(1 - S22*ZL),2))) * np.power(abs(S21),2) * (1/(1 - np.power(abs(refl_input),2)))
Gdel_dB = 10*np.log10(Gdel)

#print("Transducer power gain",GT_dB,"dB")       # mismatch at input and output
#print("Available power Gain",Gav_dB,"dB")       # mismatch at input
#print("Delivered power Gain ",Gdel_dB,"dB")     # mismatch at output

MS_dB = GT_dB - Gdel_dB
ML_dB = GT_dB - Gav_dB
MS = np.power(10,(MS_dB/10))
MS = np.power(10,(0/10))
ML = np.power(10,(ML_dB/10))
MS_max = 1/(1-np.power(abs(S11),2))
print("MS_max = ",10*np.log10(MS_max))

##gs = MS/MS_max
gs = MS*(1 - np.power(abs(S11),2))
cs = gs*np.conj(S11)/(1-np.power(abs(S11),2)*(1-gs))
rs = ((1-np.power(abs(S11),2))*np.sqrt(1 - gs))/(1-np.power(abs(S11),2)*(1-gs))
print("The centre (absolute value) of the constant gain circle in the source plane",abs(cs))
print("The centre (angle in degrees) of the constant gain circle in the source plane",np.rad2deg(cmath.phase(cs)))
print(rs)

ML_max = 1/(1-np.power(abs(S22),2))
print("ML_max = ",10*np.log10(ML_max))

gemma_L = np.conj(S22)
print("Gemma_L at the maximum value of ML",abs(gemma_L))

ML = np.power(10,(-2/10))
gl = ML*(1 - np.power(abs(S22),2))
cl = gl*np.conj(S22)/(1-np.power(abs(S22),2)*(1-gl))
rl = ((1-np.power(abs(S22),2))*np.sqrt(1 - gl))/(1-np.power(abs(S22),2)*(1-gl))
print("constant gain circle in the load plane ",abs(cl))
print("The center (angle in degrees) of the constant gain circle in the load plane",np.rad2deg(cmath.phase(cl)))
print(rl)


Gemma_opt = complex(0.27 * np.cos(np.deg2rad(75)),0.27 * np.sin(np.deg2rad(75)))
ZL_opt = 50 * (1 + Gemma_opt) / (1-Gemma_opt)
print("Determine the imaginary part of ZL_opt",ZL_opt)
L = np.imag(ZL_opt) / (2*np.pi*5.8*np.power(10,9))
print("L = ", L*np.power(10,9),"nH")
C = 1 / (2*np.pi*5.8*np.power(10,9)* (np.imag(ZL_opt)))
print("C = ", C*np.power(10,12),"pF","XXX")


## num2
Mag_S11 = 0.869
Ang_S11 = -159
Mag_S12 = 0
Ang_S12 = 0
Mag_S21 = 4.25
Ang_S21 = 61
Mag_S22 = 0.507
Ang_S22 = -117

## Convert S parameters to complex numbers
S11 = complex(Mag_S11 * np.cos(np.deg2rad(Ang_S11)),Mag_S11 * np.sin(np.deg2rad(Ang_S11)))
S21 = complex(Mag_S21 * np.cos(np.deg2rad(Ang_S21)),Mag_S21 * np.sin(np.deg2rad(Ang_S21)))
S12 = complex(Mag_S12 * np.cos(np.deg2rad(Ang_S12)),Mag_S12 * np.sin(np.deg2rad(Ang_S12)))
S22 = complex(Mag_S22 * np.cos(np.deg2rad(Ang_S22)),Mag_S22 * np.sin(np.deg2rad(Ang_S22)))
print("")
print("New")
ML_max = 1/(1-np.power(abs(S22),2))
print("ML_max = ",10*np.log10(ML_max))

ML = np.power(10,(0.5/10))
gl = ML*(1 - np.power(abs(S22),2))
cl = gl*np.conj(S22)/(1-np.power(abs(S22),2)*(1-gl))
rl = ((1-np.power(abs(S22),2))*np.sqrt(1 - gl))/(1-np.power(abs(S22),2)*(1-gl))
print("constant gain circle in the load plane ",abs(cl))
print(np.rad2deg(cmath.phase(cl)))
print(rl)


Gemma_opt = complex(0.125 * np.cos(np.deg2rad(117)),0.125 * np.sin(np.deg2rad(117)))
ZL_opt = 50 * (1 + Gemma_opt) / (1-Gemma_opt)
print("Determine the imaginary part of ZL_opt",ZL_opt)
L = np.imag(ZL_opt) / (2*np.pi*1.9*np.power(10,9))
print("L = ", L*np.power(10,9),"nH","XXX")
C = 1 / (2*np.pi*5.8*np.power(10,9)* (np.imag(ZL_opt)))
print("C = ", C*np.power(10,12),"pF","XXX")

#Length=91.8 