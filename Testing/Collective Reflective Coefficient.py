import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

ro = 1000

""" Pipe 1 """
Area1 = 100
a1 = 1000
Z1 = ro*a1/Area1

""" Pipe 2 """
Area2 = 50
a2 = 1000
Z2 = ro*a2/Area2

df = 0.001
maxf = 1000

print("Z1 =", Z1)

print("Z2 =", Z2)




frequencyRange = np.arange(0,maxf,df)
SaveMatrix = np.zeros((3,len(frequencyRange)),dtype=complex)

t = 1
x1 = 0
x2 = 1000
L= x2-x1

for f in tqdm(range(0,len(frequencyRange))):
    omega = 2*np.pi*frequencyRange[f]
    k1 = omega/a1
    k2 = omega/a2
    j=1j
    Pin = 1
    """
    P[0] = Prf1
    P[1] = Ptr1
    P[2] = Prf2
    P[3] = Ptr2
    """
    A = [[np.exp(j*(-k1*x1-omega*t)), -np.exp(j*(k2*x1-omega*t)), -np.exp(j*(-k2*x1-omega*t)), 0],
         [np.exp(j*(-k1*x1-omega*t))/Z1, np.exp(j*(k2*x1-omega*t))/Z2, -np.exp(j*(-k2*x1-omega*t))/Z2, 0],
         [0,np.exp(j*(k2*x2-omega*t)), np.exp(j*(-k2*x2 - omega*t)), -np.exp(j*(k1*(x2-L)-omega*t))],
         [0,np.exp(j*(k2*x2-omega*t))/Z2, -np.exp(j*(-k2*x2 - omega*t))/Z2, -np.exp(j*(k1*(x2-L)-omega*t))/Z1]]

    B = [[-Pin*np.exp(j*(k1*x1-omega*t))],
         [(Pin*np.exp(j*((k1*x1-omega*t))))/Z1],
         [0],
         [0]]

    P = np.linalg.solve(A,B)
    Prf1 = P[0][0]
    Ptr2 = P[3][0]
    R = Prf1 / 1
    T = Ptr2/1
    SaveMatrix[0][f] = frequencyRange[f]
    SaveMatrix[1][f] = R
    SaveMatrix[2][f] = T

f_seq = SaveMatrix[0][:]
R_seq = SaveMatrix[1][:]
T_seq = SaveMatrix[2][:]
time = np.arange(0, (1 / df), 1 / maxf)

fft_R_seq = np.real(np.fft.fft(R_seq, (len(f_seq))))
plt.plot(time, fft_R_seq)
for i in fft_R_seq:
    if abs(i) > 0.1:
        print(i)

plt.show()
