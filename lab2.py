import random
import matplotlib.pyplot as plt
import numpy as np
import math
from prettytable import PrettyTable

a = -2
b = 2


def Y(x):
    return x*x


def f(y):
    return 1/(4 * y**0.5)


def F(y):
    return y**0.5/2


def number_of_values_on_interval(M, variation, A_intervals, B_intervals):
    v_interval = [0] * M
    for i in range(0, M):
        for value in variation:
            if A_intervals[i] < value < B_intervals[i]:
                v_interval[i] += 1
            if value == B_intervals[i] and i != M - 1:
                v_interval[i] += 0.5
                v_interval[i + 1] += 0.5
    v_interval[0] += 1
    v_interval[-1] += 1
    return v_interval


def find_borders():
    ay = Y(a)
    by = Y(b)
    xes = np.arange(a, b + 0.001, 0.001)
    for x in xes:
        ay = min(ay, Y(x))
        by = max(by, Y(x))
    ay = round(ay, 2)
    by = round(by, 2)
    return ay, by


def interval_method(variation, M, n):       # равноинтервальный метод
    delta = (variation[-1] - variation[0])/M
    A_intervals = [variation[0] + (i-1)*delta for i in range(1, M+1)]
    B_intervals = [variation[0] + i*delta for i in range(1, M+1)]
    v_interval = number_of_values_on_interval(M, variation, A_intervals, B_intervals)
    f_interv = [v_interval[i]/(n*delta) for i in range(M)]  # эмпирическая функция плотности распределения
    return A_intervals, B_intervals, f_interv


def possibility_method(variation, M, n):         # равновероятностный метод
    mi = n // M
    B_pos = [(variation[i * mi] + variation[i * mi + 1]) / 2 for i in range(1, M)] + [variation[-1]]
    A_pos = [variation[0]] + B_pos[:-1]
    delta_pos = [B_pos[i] - A_pos[i] for i in range(M)]
    v_pos = number_of_values_on_interval(M, variation, A_pos, B_pos)
    f_pos = [v_pos[i] / (n * delta_pos[i]) for i in range(M)]
    return A_pos, B_pos, f_pos


def histogram_n_polygon(M, func, A, B):
    plt.title("Гистограмма\n и полигон\n распределения")
    max_height = [0] * (M + 1)
    for i in range(M):
        max_height[i] = max(func[i-1], func[i])
    max_height[0] = func[0]
    max_height[-1] = func[-1]
    plt.hlines(func, A, B, color='red')
    plt.axhline(0, color='red')
    A.append(B[-1])
    plt.vlines(A, [0], max_height, color='red')
    plt.xlim(A[0] - 0.5, B[-1] + 0.5)
    middle = [(A[i] + A[i + 1])/2 for i in range(M)]
    plt.plot(middle, func, color='blue')


def theor_density(ay, by, A, B, func):
    plt.title("Теоретическая\n плотность\n  распределения")
    xes = np.arange(ay, by, 0.001)
    ys = [f(x) for x in xes if x != 0]
    xes = [x for x in xes if x != 0]
    plt.ylim(-0.1, max(func) + 0.2)
    plt.xlim(A[0] - 0.5, B[-1] + 0.5)
    plt.hlines(0, A[0] - 0.5, 0, color='blue')
    plt.hlines(0, 4, B[-1] + 0.5, color='blue')
    plt.plot(xes, ys, color='blue')


def emp_function(func, A, B, M):
    plt.title("Эмпирическая\n функция\n  распределения")
    F_emp = [0] * M
    F_emp[0] = func[0]*(B[0] - A[0])
    for i in range(1, M):
        F_emp[i] += F_emp[i-1] + func[i] * (B[i] - A[i])
    plt.hlines(F_emp, A, B, lw=2)
    print("Empirical function")
    f_emp_table = PrettyTable()
    f_emp_table.field_names = ["Interval", "f emp"]
    for i in range(len(func)):
        f_emp_table.add_row(["%.4f" % A[i] + " : " + "%.4f" % B[i], func[i]])
    print(f_emp_table)


def graphics(M, func, A, B):
    plt.subplot(1, 3, 1)
    histogram_n_polygon(M, func, A, B)
    plt.grid()
    plt.subplot(1, 3, 2)
    ay, by = find_borders()
    theor_density(ay, by, A, B, func)
    plt.grid()
    plt.subplot(2, 3, 3)
    emp_function(func, A, B, M)
    plt.grid()
    plt.show()


if __name__ == "__main__":
    n = int(input("Enter n.\n"))
    eps = [np.random.uniform(0, 1) for i in range(n)]
    x = [(ei*(b-a) + a) for ei in eps]
    y = [Y(xi) for xi in x]
    variation = sorted(y)
    if n <= 100:
        M = int(n**0.5)
    else:
        M = int(random.uniform(2, 4) * math.log10(n))
    print_variation = [round(v, 6) for v in variation]
    print(print_variation)
    A_intervals, B_intervals, f_interv = interval_method(variation, M, n)
    A_pos, B_pos, f_pos = possibility_method(variation, M, n)
    print("Равноинтервальный метод: \n")
    graphics(M, f_interv, A_intervals, B_intervals)
    print("Равновероятностный метод: \n")
    graphics(M, f_pos, A_pos, B_pos)


