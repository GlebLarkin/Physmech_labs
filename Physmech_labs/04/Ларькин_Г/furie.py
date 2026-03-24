import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# NewFile13_32hz.csv

data = []
with open("./NewFile13_32hz.csv", 'r') as f:
    f.readline()
    f.readline()
    lines = f.readlines()
    for line in lines:
        data.append([float(i) for i in line[:-2].strip().split(',')])
data = np.array(data)


def furie(data):
    N = len(data)
    dt = 2e-3
    samplerate = 1 / dt

    ff = 2 / N * np.fft.fft(data[:, 1])[:N // 2]
    ffOUT = 2 / N * np.fft.fft(data[:, 1])
    ffFreq = np.fft.fftfreq(N, dt)[:N // 2]

    plt.plot(ffFreq, np.abs(ff), 'o')
    plt.axvline(x=3.33, color='red', linestyle='--')


    points_to_label = [8, 16]

    for idx in points_to_label:
        plt.annotate(f'{ffFreq[idx]:.2f} Гц',
                     xy=(ffFreq[idx], np.abs(ff)[idx]),
                     xytext=(ffFreq[idx] + 0.2, np.abs(ff)[idx]),
                     fontsize=9,
                     ha='left', va='center',
                     bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7))


    plt.scatter(ffFreq[points_to_label], np.abs(ff)[points_to_label],
                c='red', s=30, zorder=5)

    plt.minorticks_on()
    plt.xlim(0, 8)


    plt.xlabel('Частота (Гц)', fontsize=12)
    plt.ylabel('Амплитуда (В)', fontsize=12)
    plt.title('Спектр сигнала', fontsize=14)

    plt.grid(which='minor', color="gray", lw=0.5, ls='--', alpha=0.3)
    plt.grid(which='major', color="black", lw=0.5, alpha=0.5)
    plt.show()


furie(data)