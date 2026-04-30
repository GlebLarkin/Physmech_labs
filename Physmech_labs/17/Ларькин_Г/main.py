import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

T = 293
P = 739
M = 0.029
R = 8.31
CALIBER = 6e-2
DIST_BETW_PIPES = 1e-2

def Torr2Pa(torr):
    return torr * 133.322

def Pa2Speed(pa):
    return np.sqrt((2 * pa * R * T) / (Torr2Pa(P) * M))

#ro = p * mu /(R T)
df = pd.read_excel('data.ods', engine='odf')
data = df.to_numpy(dtype=int)

n_calibers, n_cols = data.shape
atm_per_row = (data[:, 0] + data[:, -1]) / 2.0
center_idx = n_cols // 2
coords = (np.arange(n_cols) - center_idx) * DIST_BETW_PIPES * 100  # см

caliber_labels = [0, 2, 4, 6, 8, 10, 12.5, 15]
x_over_D = np.array(caliber_labels)  # расстояние в калибрах

cmap = mpl.colormaps['plasma']
colors = [cmap(i / (n_calibers - 1)) for i in range(n_calibers)]

speeds_list = []
u_max_list = []
r_12_list = []

for i in range(n_calibers):
    #dyn_pa = Torr2Pa(atm_per_row[i] - data[i, :])
    dyn_pa = (atm_per_row[i] - data[i, :]) * 10
    dyn_pa = np.clip(dyn_pa, 0, None)
    speed = Pa2Speed(dyn_pa)
    speeds_list.append(speed)

    u_max = speed[center_idx]
    u_max_list.append(u_max)
    half = u_max / 2.0

    right_coords = coords[center_idx:]
    right_speed = speed[center_idx:]
    r_12 = None
    for j in range(len(right_speed) - 1):
        if right_speed[j] >= half >= right_speed[j+1]:
            x0, x1 = right_coords[j], right_coords[j+1]
            y0, y1 = right_speed[j], right_speed[j+1]
            r_12 = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            break
    r_12_list.append(r_12)

u_max_arr = np.array(u_max_list)
r_12_arr  = np.array(r_12_list)
u_max_0   = u_max_arr[0]  # скорость на оси в первом сечении

fig1, ax1 = plt.subplots(figsize=(10, 6))
for i in range(n_calibers):
    ax1.plot(coords, speeds_list[i], marker='o', markersize=4,
             color=colors[i], label=f'd = {caliber_labels[i]}D')
    ax1.axhline(u_max_list[i]/2, color=colors[i], ls='--', lw=0.8, alpha=0.6)
    if r_12_arr[i] is not None:
        ax1.axvline( r_12_arr[i], color=colors[i], ls=':', lw=0.8, alpha=0.6)
        ax1.axvline(-r_12_arr[i], color=colors[i], ls=':', lw=0.8, alpha=0.6)
ax1.set_xlabel('Координата x, см', fontsize=13)
ax1.set_ylabel('Скорость, м/с', fontsize=13)
ax1.set_title('Профили скорости (пунктир — U_max/2, точки — ±r½)', fontsize=14)
ax1.legend(title='Расстояние до сопла', fontsize=9)
ax1.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('velocity_profile.png', dpi=150)
plt.show()

fig2, ax2 = plt.subplots(figsize=(10, 6))
for i in range(n_calibers):
    if r_12_arr[i] is None:
        continue
    r_norm = coords / r_12_arr[i]
    u_norm = speeds_list[i] / u_max_arr[i]
    ax2.plot(r_norm, u_norm, marker='o', markersize=4,
             color=colors[i], label=f'd = {caliber_labels[i]}D')
ax2.axhline(0.5, color='gray', ls='--', lw=1, label='U/U_max = 0.5')
ax2.axvline( 1.0, color='gray', ls=':', lw=1)
ax2.axvline(-1.0, color='gray', ls=':', lw=1)
ax2.set_xlabel('r / r½', fontsize=13)
ax2.set_ylabel('U / U_max', fontsize=13)
ax2.set_title('Нормированные профили скоростей', fontsize=14)
ax2.legend(title='Расстояние до сопла', fontsize=9)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('velocity_normalized.png', dpi=150)
plt.show()

fig3, axes = plt.subplots(1, 2, figsize=(13, 5))

ax3 = axes[0]
ax3.scatter(x_over_D, 1 / (u_max_arr / u_max_0), color='steelblue')
ax3.set_xlabel('(x / D)^-1', fontsize=13)
ax3.set_ylabel('U_max(x) / U_max(0)', fontsize=13)
ax3.set_title('Падение скорости на оси', fontsize=14)
ax3.grid(True, alpha=0.3)

# Линеаризация: теория U_max ~ 1/x  =>  1/U_max ~ x
# Фитируем только основной участок (пропускаем x=0, берём x>0)
mask = x_over_D > 0
x_fit = x_over_D[mask]
y_fit = 1.0 / (u_max_arr[mask] / u_max_0)  # U_max(0)/U_max(x)

ax3b = ax3.twinx()
#ax3b.plot(x_over_D[mask], y_fit, 's--', color='tomato', lw=1.5, ms=6, label='U_max(0)/U_max')
coeffs = np.polyfit(x_fit, y_fit, 1)
x_line = np.linspace(x_fit[0], x_fit[-1], 100)
#ax3b.plot(x_line, np.polyval(coeffs, x_line), '-', color='tomato', alpha=0.5,
          #label=f'Линейный фит: {coeffs[0]:.3f}·x/D + {coeffs[1]:.3f}')
#ax3b.set_ylabel('U_max(0) / U_max(x)  (правая ось)', fontsize=11, color='tomato')
#ax3b.tick_params(axis='y', labelcolor='tomato')
#ax3b.legend(fontsize=9, loc='upper left')

ax4 = axes[1]
ax4.scatter(x_over_D, r_12_arr, color='seagreen')

# Линейный фит (теория: r½ = C1 * x * D, т.е. r½/D ~ x/D)
# r_12 в см, D = 6 см
D_cm = CALIBER * 100
r12_over_D = r_12_arr / D_cm
coeffs2 = np.polyfit(x_over_D[mask], r12_over_D[mask], 1)
x_line2 = np.linspace(x_over_D[mask][0], x_over_D[mask][-1], 100)
ax4.plot(x_line2, np.polyval(coeffs2, x_line2) * D_cm, '--', color='seagreen',
         alpha=0.6, label=f'Аппроксимация: C₁ = {coeffs2[0]:.3f}')

ax4.set_xlabel('x / D', fontsize=13)
ax4.set_ylabel('r½, см', fontsize=13)
ax4.set_title('Рост полуширины струи', fontsize=14)
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3)

plt.suptitle('Осевые характеристики струи', fontsize=15, y=1.01)
plt.tight_layout()
plt.savefig('characteristics.png', dpi=150, bbox_inches='tight')
plt.show()

print("C1 (наклон r½/D от x/D):", coeffs2[0])
print("Коэффициенты фита 1/U_max:", coeffs)
