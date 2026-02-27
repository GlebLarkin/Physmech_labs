import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4
mu = 28.9e-3
R = 8.31


def pi(M):
    return (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M * M)) ** (-gamma / (gamma - 1))


def tau(M):
    return 1 / (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M * M))


def q(M):
    return M * (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M * M)) ** (-(gamma + 1) / (2 * (gamma - 1)))


def M6_(q):
    M1 = 1
    M2 = 6

    def F(M):
        q_tmp = M * (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M * M)) ** (-(gamma + 1) / (2 * (gamma - 1)))
        return q_tmp - q

    for i in range(1, 100):
        M = (M1 + M2) / 2
        # print(M)

        if (F(M2) * F(M) < 0):
            M1 = M
        else:
            M2 = M

    return M


def M3_(q):
    M1 = 0
    M2 = 1

    def F(M):
        q_tmp = M * (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M * M)) ** (-(gamma + 1) / (2 * (gamma - 1)))
        return q_tmp - q

    for i in range(1, 100):
        M = (M1 + M2) / 2
        # print(M)

        if (F(M2) * F(M) < 0):
            M1 = M
        else:
            M2 = M

    return M


def fM(M):
    return np.sqrt((1 + (gamma - 1) / 2 * M * M) / (gamma * M * M - (gamma - 1) / 2))
    # return np.sqrt(1 / (((gamma + 1) / 2) **2 * M*M / (1 + (gamma - 1) / 2 * M*M) + (gamma - 1) / 2))


def fP(M):
    return 2 * gamma / (gamma + 1) * M * M - (gamma - 1) / (gamma + 1)


def fT(M):
    return (2 * gamma * M * M - gamma + 1) * ((gamma - 1) * M * M + 2) / ((gamma + 1) * M) ** 2


def u2a1_(p2p1):
    return (p2p1 - 1) / np.sqrt((gamma * (gamma - 1) / 2) * (1 + (gamma + 1) / (gamma - 1) * p2p1))


def p2p1_(u2a1):
    p2p1 = 1
    for _ in range(100):
        p2p1 = u2a1 * np.sqrt((gamma * (gamma - 1) / 2) * (1 + (gamma + 1) / (gamma - 1) * p2p1)) + 1
        # print(p2p1)
    return p2p1


def p2p1__(M_sh):
    return (M_sh ** 2 - (gamma - 1) / 2 / gamma) * gamma * 2 / (gamma + 1)


def main():
    with open('results.csv', 'w', newline='', encoding='utf-8') as csvfile:

        x1 = []
        y1 = []
        for p5p4 in np.linspace(1, 0.279, 100):
            F4 = 1
            F1 = 5
            a5a4 = p5p4 ** ((gamma - 1) / 2 / gamma)
            M5 = (1 / a5a4 - 1) * 2 / (gamma - 1)
            q3 = F4 / F1 * q(M5)
            M3 = M3_(q3)
            p3p5 = pi(M3) / pi(M5)
            f = p3p5 * p5p4

            a4a1 = 1
            a3a5 = np.sqrt(tau(M3) / tau(M5))
            u2u1 = M3 * a4a1 * a3a5 * a5a4
            p2p1 = p2p1_(u2u1)
            y1.append(p2p1)
            p4p1 = p2p1 / f
            x1.append(p4p1)

        x2 = []
        y2 = []
        for Fs in np.linspace(1, 5, 100):
            # Fs = np.linspace(3, 5, 1000)
            M5 = 1
            F4 = 1  # критической сечение
            F1 = 5
            q6 = F4 / Fs * q(M5)
            M6 = M6_(q6)
            M7 = fM(M6)
            q3 = Fs / F1 * q(M7)
            M3 = M3_(q3)
            # print(M6)

            p3p7 = pi(M3) / pi(M7)
            p7p6 = fP(M6)
            p6p5 = pi(M6) / pi(M5)
            p5p4 = (2 / (gamma + 1)) ** (2 * gamma / (gamma - 1))

            f = p3p7 * p7p6 * p6p5 * p5p4

            a4a1 = 1
            a3a7 = np.sqrt(tau(M3) / tau(M7))
            a7a6 = np.sqrt(fT(M6))
            a6a5 = np.sqrt(tau(M6) / tau(M5))
            a5a4 = 2 / (gamma + 1)
            u2a1 = M3 * a4a1 * a3a7 * a7a6 * a6a5 * a5a4
            p2p1 = p2p1_(u2a1)
            y2.append(p2p1)
            p4p1 = p2p1 / f
            x2.append(p4p1)

    ############################
    T = 297
    t1 = np.array([198, 191, 186, 182, 178, 176, 170]) * 1e-6
    t2 = np.array([330, 320, 312, 305, 297, 293, 284]) * 1e-6
    p4_au = np.array([125, 125, 125, 167, 125, 153, 125])
    p4 = 9.81e4 + p4_au * 6 * 9.81e4 / 250
    p1 = np.array([752, 500, 375, 375, 250, 250, 188]) * 9.81e4 / 752
    p4p1_exp = p4 / p1
    p2p1_exp1 = []
    l = 5e-2
    u_sh = l / (t2 - t1)
    a1 = np.sqrt(gamma * R * T / mu)
    print(a1)
    for u_sh_ in u_sh:
        M_sh = u_sh_ / a1
        p2p1_exp1.append(p2p1__(M_sh))

    dV = np.array([91.6, 89.2, 83.6, 95.6, 74.8, 83.2, 69.2])
    k = 3.6e-3
    p2p1_exp2 = 1 + dV / k / p1

    print(p2p1_exp2)

    plt.plot(x1, y1, c='b', label='теория без скачка')
    plt.plot(x2, y2, '--', c='b', label='теория со скачком')
    plt.scatter(p4p1_exp, np.array(p2p1_exp1), c='purple', label='эксперимент(по значению скорости)')
    plt.scatter(p4p1_exp, np.array(p2p1_exp2), c='g', label='эксперимент(по значению напряжения)')

    plt.xlim([1, 17])
    plt.ylim([1, 2])

    plt.xlabel('p4/p1')
    plt.ylabel('p2/p1')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    plt.legend()
    plt.show()
    return


if __name__ == '__main__':
    main()
