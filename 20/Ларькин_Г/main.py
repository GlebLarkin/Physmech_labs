''' 
Лабораторная работа 20
*Гидродинамическая устойчивость вращательного движения течения Куэтта*
Авторы - Глеб Ларькин, Валентина Копышева, Дай Сюй, Себастиан Хулио
11.02.26


Вращаются два соосных цилиндра, между ними жидкость с алюминиевой пудрой.
Рассматриваются три режима:
- вращение с одинакоывми частотами (ничего не измеряется), 
- вращение с малой разностью угловых скоростей (ничего не измеряется),
- вращение с большой разностью угловых скоростей (течение Куэтта).


При разных частотах вращения устанавливаются
разные течения:

- В первом режиме линейное распределение скоростей.

- Во втором режиме распределение скоростей имеет вид 
u(r) = A r + B / r (см лекции Жмура/лабник).

В первых двух режимах, несмотря на градиент давлений, вязкость уравнивает силы
давлений, и течение устойчивое.

- В третьем режиме градиент давлений становится больше вязких сил, 
и начинается также и перпендикулярное течение. 
Итоговое течение - течение Куэтта - сумма двух движений по окружностей
во взаимноперпендикулярных направлениях. 

Такое течение - неустойчивое, но стационарное. 
Ожидаем, что все точки будут в области неустойчивости,
но чем меньше вихрей, тем ближе мы к области устойчивости.

Измеряется число вихрей в зависимости от частоты вращения.

!Важно: направленность течения в вихрях чередуется, поэтому некоторые полосы 
(например, четные) видно лучше других (например, нечетных).

Жидкость прилипает к стенкам, поэтому число Рейндольдса течения
жидкости рядом с цилиндром определяется скоростью цилиндра:
    Re_i = 2 pi f r_i^2 / nu

Характерный размер вихря = разность радиусов - размер пограничных слоев.

В работе также исследуется вращение цилиндров в противоположные стороны:
в этом режиме скорость в какой-то момент меняет знак => 
существует поверхность нулевой скорости. Эта поверхность - окружность, 
так как иначе была бы несимметрия, и возникала бы неустойчивость, 
такое течение установиться не могло.
Тогда можно заменить исходную задачу на две независимые, 
где одна стенка цилиндра вращается, 
а вторая покоится => будет два набора вихрей.

Система обладает гистерезисом, результаты плохо повторяемы.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp 
import pandas as pd

def get_Re(f, r):
    '''
    Возвращает число Рейнольдса при заданной циклической частоте и радиусе
    '''
    return 2 * np.pi * f * r**2 / consts.nu

def get_Re_error(f, r, delta_f=1.0):
    """
    Возвращает абсолютную погрешность числа Рейнольдса.
    """
    Re = get_Re(f, r)
    
    dRe_df = 2 * np.pi * r**2 / consts.nu
    term_f = (dRe_df * delta_f)**2
    
    delta_nu = consts.nu_error
    term_nu = (Re / consts.nu * delta_nu)**2
    
    delta_Re = np.sqrt(term_f + term_nu)
    
    return delta_Re

#############################################

class consts:
    r1 = 3 * 10**(-2)   # м (внутренний цилиндр)
    r2 = 3.5 * 10**(-2) # м (внешний цилиндр)
    
    nu = 2 * 10**(-5) # м^2/с - кинематическая вязкость
    nu_error = nu * 0.05 #  м^2/с - погрешность для кинематической вязкости
    
    Re_df = pd.DataFrame(
    {
        "Re2": [-4000, -3500, -3000, -2500, -2000, -1500, -1000, -500, 0,
                 0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000],
        "Re1": [2055, 1935, 1760, 1575, 1390, 1185, 975, 785, 632,
                 632, 833, 1206, 1645, 2111, 2590, 3076, 3565, 4057],
        "Вращение": ["Встречное"] * 9 + ["Параллельное"] * 9
    }) # Граница устойчивости - табличные значения

class data_exp1:
    ### Первая серия экспериментов - внешний цилиндр не вращается
    
    # количество вихрей
    harmonics = np.array([1, 3, 10, 14, 20, 24, 34]) #! Данные измерены очень криво
    f1 = np.array([3.4, 5.7, 8.0, 8.1, 8.3, 8.6, 8.8]) # Гц (Показания частометра / 100)
    f2 = np.array([0] * f1.size) # Гц (Цилиндр не вращался)
    
    df = pd.DataFrame({'Number of harmonics': harmonics, 'f1': f1, 'f2': f2})
    
    df['Re1'] = df['f1'].apply(lambda x: get_Re(x, consts.r1))
    df['Re2'] = df['f2'].apply(lambda x: get_Re(x, consts.r2))
    
    df['Re1_err'] = [get_Re_error(f, consts.r1) for f in df['f1']]
    df['Re2_err'] = [get_Re_error(f, consts.r2) for f in df['f2']]

class data_exp2:
    ### Вторая серия экспериментов - внешний цилиндр вращается
    ### с частотой ПРИМЕРНО 1 Гц 
    ### (у нас сломался частотометр, мы смотрли таймером)

    # количество вихрей
    harmonics = np.array([3, 12, 18, 36])
    # Смотрим вращение относительно внешнего цилиндра 
    f1 = np.array([8.8, 9.5, 9.9, 10]) # Гц (Показания частометра / 100)
    f2 = np.array([1] * f1.size) # Гц !с ОГРОМНОЙ погрешностью
    
    df = pd.DataFrame({'Number of harmonics': harmonics, 'f1': f1, 'f2': f2})
    df['Re1'] = df['f1'].apply(lambda x: get_Re(x, consts.r1))
    df['Re2'] = df['f2'].apply(lambda x: get_Re(x, consts.r2))
    
    df['Re1_err'] = [get_Re_error(f, consts.r1) for f in df['f1']]
    df['Re2_err'] = [get_Re_error(f, consts.r2) for f in df['f2']]
    
class data_exp3:
    ### Третья серия экспериментов - вращение цилиндров в разные стороны
    f1 = 4.25 # Гц  (Показания частотометра)
    f2 = -0.4 # Гц (Замеры с помощью секундомера)
    # При этих частотах Установился режим в двумя вихрями
    Re1 = get_Re(f1, consts.r1)
    Re2 = get_Re(f2, consts.r2)
    
    Re1_err = get_Re_error(f1, consts.r1)
    Re2_err = get_Re_error(f2, consts.r2)

print(consts.Re_df)
print(data_exp1.df)
print(data_exp2.df)
print(data_exp3.Re1, data_exp3.Re2)

# Re1(Re2)
df_consts = consts.Re_df

plt.plot(df_consts["Re2"], df_consts["Re1"], 
         'r-', linewidth=2, label='Граница устойчивости')

plt.errorbar(data_exp1.df["Re2"], data_exp1.df["Re1"], 
             xerr=data_exp1.df['Re2_err'], yerr=data_exp1.df['Re1_err'],
             fmt='o', color='green', markersize=8, capsize=3, label='Эксп. 1 (f2=0)')

plt.errorbar(data_exp2.df["Re2"], data_exp2.df["Re1"], 
             xerr=data_exp2.df['Re2_err'], yerr=data_exp2.df['Re1_err'],
             fmt='s', color='orange', markersize=8, capsize=3, label='Эксп. 2 (f2≈1 Гц)')

plt.errorbar(data_exp3.Re2, data_exp3.Re1, 
             xerr=data_exp3.Re2_err, yerr=data_exp3.Re1_err,
             fmt='D', color='purple', markersize=8, capsize=3, label='Эксп. 3 (встречное)')

plt.xlabel(r'$Re_2$', fontsize=14)
plt.ylabel(r'$Re_1$', fontsize=14)
plt.title('Зависимость $Re_1$ от $Re_2$ (граница устойчивости и эксперименты)', fontsize=14)
plt.minorticks_on()
plt.grid(True, which='major', linestyle='--', alpha=0.7)
plt.grid(True, which='minor', linestyle=':', alpha=0.4)
plt.legend(fontsize=12)
plt.tight_layout()

# f(harmonics)
plt.figure()
plt.scatter(data_exp1.df['f1'], data_exp1.df['Number of harmonics'],
            color='green', s=80, marker='o')
plt.ylabel('Количество вихрей', fontsize=12)
plt.xlabel('f1, Гц', fontsize=12)
plt.title('Зависимость числа вихрей от частоты\n(внешний цилиндр неподвижен)', fontsize=14)
plt.grid(True, which='major', linestyle='--', alpha=0.7)
plt.grid(True, which='minor', linestyle=':', alpha=0.4)
plt.tight_layout()

plt.figure()
delta_f = data_exp2.df['f1'] - data_exp2.df['f2']
plt.scatter(delta_f, data_exp2.df['Number of harmonics'],
            color='orange', s=80, marker='s')
plt.ylabel('Количество вихрей', fontsize=12)
plt.xlabel('f1 - f2, Гц', fontsize=12)
plt.title('Зависимость числа вихрей от разности частот\n(внешний цилиндр вращается с f2 ≈ 1 Гц)', fontsize=14)
plt.grid(True, which='major', linestyle='--', alpha=0.7)
plt.grid(True, which='minor', linestyle=':', alpha=0.4)
plt.tight_layout()

plt.show()
