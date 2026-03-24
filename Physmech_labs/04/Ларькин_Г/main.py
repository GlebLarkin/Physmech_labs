import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sp
import pandas as pd

class consts:
    h     = 3.5  * 1e-2  # м — глубина жидкости в кювете
    dh    = 0.5  * 1e-2  # м — погрешность h
    g     = 9.8          # м/с²
    L     = 1.44         # м — длина кюветы
    dL    = 3.0  * 1e-2  # м — погрешность L

class data_exp1:
    ### Зависимость omega от k

    # ── Теория ──────────────────────────────────────────────────────────────
    n_theor      = np.array(range(1, 30))
    kappa_theor  = np.pi * n_theor / consts.L
    freq_theor   = np.sqrt(kappa_theor * consts.g * np.tanh(kappa_theor * consts.h)) / (2 * np.pi)
    omega_theor  = 2 * np.pi * freq_theor

    # Погрешность теоретической кривой по h
    # d(omega)/d(h) через дисперсионное соотношение:
    # omega = sqrt(g*k*tanh(k*h))  =>
    # d(omega)/d(h) = g*k^2 / (2*omega * cosh^2(k*h))
    domega_theor_dh = (consts.g * kappa_theor**2
                       / (2 * omega_theor * np.cosh(kappa_theor * consts.h)**2))
    domega_theor    = domega_theor_dh * consts.dh   # абс. погрешность omega_theor

    # ── Эксперимент ─────────────────────────────────────────────────────────
    # Гц — выставляемые частоты генератора (теоретически предсказанные)
    freq_for_exp = np.array([1.85, 2.0, 2.15, 2.29, 2.43, 2.56,
                             2.68, 2.80, 2.91, 3.02, 3.13])
    # м — измеренные длины волн
    lambdas_exp  = np.array([24, 18, 17, 13, 13, 15,
                             14, 14, 11, 12, 11]) * 1e-2
    dlambda      = 3e-2  # м — погрешность lambda

    omega_exp  = 2 * np.pi * freq_for_exp
    kappa_exp  = 2 * np.pi / lambdas_exp

    # Погрешность kappa: kappa = 2*pi/lambda  =>  d(kappa) = kappa * dlambda / lambda
    dkappa_exp = kappa_exp * dlambda / lambdas_exp


# ── График ──────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(8, 5.5))
fig.patch.set_facecolor("#FAFAF7")
ax.set_facecolor("#FAFAF7")

# Теоретическая кривая с полосой погрешности
ax.plot(data_exp1.kappa_theor, data_exp1.omega_theor,
        color="#1a3a5c", lw=2, label=r"Теория: $\omega = \sqrt{g\kappa\tanh(\kappa h)}$", zorder=3)
ax.fill_between(data_exp1.kappa_theor,
                data_exp1.omega_theor - data_exp1.domega_theor,
                data_exp1.omega_theor + data_exp1.domega_theor,
                color="#1a3a5c", alpha=0.15, label=r"Погрешность по $h$", zorder=2)

# Экспериментальные точки с крестами погрешностей
ax.errorbar(data_exp1.kappa_exp, data_exp1.omega_exp,
            xerr=data_exp1.dkappa_exp,
            fmt="o", color="#c0392b", markersize=6,
            ecolor="#c0392b", elinewidth=1.4, capsize=4, capthick=1.4,
            label=r"Эксперимент $\pm\,\delta\kappa$", zorder=4)

# Оформление
ax.set_xlabel(r"Волновое число $\kappa$, м$^{-1}$", fontsize=13)
ax.set_ylabel(r"Циклическая частота $\omega$, рад/с", fontsize=13)
ax.set_title(r"Дисперсионное соотношение: $\omega(\kappa)$", fontsize=14, pad=12)

ax.legend(framealpha=0.6, fontsize=11, loc="upper left")
ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.5, color="#888")
ax.tick_params(labelsize=11)

fig.tight_layout()

plt.show()
