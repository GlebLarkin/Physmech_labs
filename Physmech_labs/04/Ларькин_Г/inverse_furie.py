import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = []
with open("./NewFile13_32hz.csv", 'r') as f:
    f.readline()
    f.readline()
    lines = f.readlines()
    for line in lines:
        data.append([float(i) for i in line[:-2].strip().split(',')])
data = np.array(data)


def inv_furie(data):

    N = len(data)
    dt = 2e-3 
    t = np.arange(N) * dt 
    X = 2 / N * np.fft.fft(data[:, 1])[:N // 2]

    a0 = X[0]
    a8 = X[8]
    a16 = X[16]

    print(a0, a8, a16)

    f0 = 0
    f8 = 8 / (N * dt)
    f16 = 16 / (N * dt)

    y_recon = a0 * np.exp(1j * 2 * np.pi * f0 * t) + \
              a8 * np.exp(1j * 2 * np.pi * f8 * t) + \
              a16 * np.exp(1j * 2 * np.pi * f16 * t)

    y_recon2 = a0 * np.exp(1j * 2 * np.pi * f0 * t) + \
              a8 * np.exp(1j * 2 * np.pi * f8 * t)

    plt.figure(figsize=(12, 5))
    plt.plot(t, data[:, 1], alpha=0.7, label = 'сигнал')
    plt.plot(t, y_recon.real, '--', linewidth=2, label = 'основная гармоника + вспомогательная')
    plt.plot(t, y_recon2.real, '--', c= 'r',  linewidth=2, label = 'основная гармоника')
    plt.xlabel('Время (с)')
    plt.ylabel('Амплитуда')

    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()



inv_furie(data)
