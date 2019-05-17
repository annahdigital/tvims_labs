import matplotlib.pyplot as plt
import numpy as np
from prettytable import PrettyTable

a = -2
b = 2


def Y(x):
    return x*x


def F(y):
    return y**0.5/2


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


def emp_theor_graphics(ay, by, variation):
    hi = (by - ay) / n
    h_start = [hi*i for i in range(n)]
    h_end = [h + hi for h in h_start]
    f_emp = [len([v for v in variation if v <= hi])/n for hi in h_end]
    plt.subplot(1, 2, 1)
    plt.hlines(f_emp, h_start, h_end, color='red')
    plt.xlim(ay, by)
    plt.title("Эмпирическая функция\n распределения")
    plt.subplot(1, 2, 2)
    arg_list = np.arange(ay, by, 0.001)
    func_theor = [F(x) for x in arg_list]
    plt.plot(arg_list, func_theor, color='blue')
    plt.xlim(ay, by)
    plt.title("Теоретическая функция\n распределения")
    plt.grid()
    plt.show()
    return h_start, h_end, f_emp


if __name__ == "__main__":
    n = int(input("Enter n.\n"))
    eps = [np.random.uniform(0, 1) for i in range(n)]
    x = [(ei*(b-a) + a) for ei in eps]
    y = [Y(xi) for xi in x]
    variation = sorted(y)
    ay, by = find_borders()
    start, end, f_emp = emp_theor_graphics(ay, by, variation)
    print("Empirical function")
    f_emp_table = PrettyTable()
    f_emp_table.field_names = ["Interval", "f emp"]
    for i in range(len(f_emp)):
        f_emp_table.add_row(["%.4f" % start[i] + " : " + "%.4f" % end[i], f_emp[i]])
    print(f_emp_table)



