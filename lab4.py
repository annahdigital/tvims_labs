import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.special import erf
from scipy.stats import norm


a = -2
b = 2
mo = 4/3
dispersion = 64/45
s = dispersion ** 0.5


def Y(x):
    return x*x


def generation_of_the_row(n):
    eps = [np.random.uniform(0, 1) for i in range(n)]
    x = [(ei * (b - a) + a) for ei in eps]
    y = [Y(xi) for xi in x]
    variation = sorted(y)
    return variation


def confidence_interval_for_expected_value(n, a, mo_es, s_e):
    print("Доверительный интервал для матожидания со значимостью ", 1 - a, ": ",
          "%.4f" % (mo_es - s_e*t.ppf(1 - a / 2, n - 1) / n**0.5), " ≤ mₓ ≤ ",
          "%.4f" % (mo_es + s_e*t.ppf(1 - a / 2, n - 1) / n**0.5))


def confidence_interval_for_expected_value_with_dispersion(n, a, mo_es):
    print("Доверительный интервал для матожидания с известной дисперсией со значимостью ", 1 - a, ": ",
          "%.4f" % (mo_es - s * norm.ppf(1-a/2) / n ** 0.5), " ≤ mₓ ≤ ",
          "%.4f" % (mo_es + s * norm.ppf(1-a/2) / n ** 0.5))


def confidence_intervals_for_expected_values(n, mo_es, s_e):
    a_values = [0.01, 0.02, 0.05]
    for a in a_values:
        confidence_interval_for_expected_value(n, a, mo_es, s_estimated)
    for a in a_values:
        confidence_interval_for_expected_value_with_dispersion(n, a, mo_es)


def confidence_interval_dependency(n, s_es):
    a_list = np.arange(0, 1, 0.001)
    interval_list = [2 * s_es * t.ppf(1 - a / 2, n - 1) / n ** 0.5 for a in a_list]
    interval_list_with_dispersion = [2 * s / n ** 0.5 * norm.ppf(1-a/2.0) for a in a_list]
    plt.plot(1 - a_list, interval_list, label='Без известной дисперсии')
    plt.plot(1 - a_list, interval_list_with_dispersion, label='С известной дисперсией')
    plt.title("Зависимость величины доверительного интервала\n матожидания от уровня значимости")
    plt.legend()
    plt.show()


def confidence_interval_volume_dependency():
    laplace_function = lambda x: erf(x / 2 ** 0.5)
    a = 0.01
    n_list = np.arange(10, 251, 1)
    interval_list = []
    interval_list_with_dispersion = []
    for n in n_list:
        var = generation_of_the_row(n)
        mo_es = sum(var) / n
        dispersion_es = sum([(x - mo_es) ** 2 for x in var]) / (n - 1)
        interval_list.append(2 * dispersion_es ** 0.5 * t.ppf(1 - a / 2, n - 1) / (n) ** 0.5)
        interval_list_with_dispersion.append(2 * s / n**0.5 * norm.ppf(1 - a / 2.0))
    plt.plot(n_list, interval_list, label='Без известной дисперсии')
    plt.plot(n_list, interval_list_with_dispersion, label='С известной дисперсией')
    plt.title("Зависимость величины доверительного интервала матожидания\n"
              " от объёма выборки с доверительным значением " + str(1 - a))
    plt.legend()
    plt.show()


def confidence_interval_for_dispersion(n, a, s_es):
    print("Доверительный интервал для дисперсии со значимостью ", 1 - a, ": ",
          "%.4f" % ((n - 1) * s_es**2 / chi2.isf((1 - (1 - a)) / 2, n - 1)), " ≤ Dₓ ≤ ",
          "%.4f" % ((n - 1) * s_es**2 / chi2.isf((1 + (1 - a)) / 2, n - 1)))


def confidence_interval_for_dispersion_with_estimation(n, a, var):
    dispersion_es = (sum([(x - mo)**2 for x in var]) / (n - 1))
    print("Доверительный интервал для дисперсии с известным матожиданием со значимостью ", 1 - a, ": ",
          "%.4f" % (n * dispersion_es / chi2.isf((1 - (1 - a)) / 2, n)), " ≤ Dₓ ≤ ",
          "%.4f" % (n*dispersion_es / chi2.isf((1 + (1 - a)) / 2, n)))


def confidence_interval_for_all_dispersions(n, s_es, var):
    a_values = [0.01, 0.02, 0.05]
    for a in a_values:
        confidence_interval_for_dispersion(n, a, s_es)
    for a in a_values:
        confidence_interval_for_dispersion_with_estimation(n, a, var)


def confidence_interval_dependency_for_dispersion(n, var, d_es):
    a_list = np.arange(0.001, 1, 0.001)
    interval_list = [(n - 1) * d_es / chi2.isf((1 + (1 - a)) / 2, n - 1) -
                     (n - 1) * d_es / chi2.isf((1 - (1 - a)) / 2, n - 1) for a in a_list]
    dispersion_es_with_mo = sum([(x - mo) ** 2 for x in var]) / (n - 1)
    interval_list_with_mo = [n*dispersion_es_with_mo / chi2.isf((1 + (1 - a)) / 2, n) -
                     n * dispersion_es_with_mo / chi2.isf((1 - (1 - a)) / 2, n) for a in a_list]
    plt.plot(1 - a_list, interval_list, label='Без известного матожидания')
    plt.plot(1 - a_list, interval_list_with_mo, label='С известным МО')
    plt.title("Зависимость величины доверительного интервала\n дисперсии от уровня значимости")
    plt.legend()
    plt.show()


def confidence_interval_volume_dependency_for_dispersion():
    a = 0.01
    n_list = np.arange(15, 151, 1)
    interval_list = []
    interval_list_with_mo = []
    for n in n_list:
        var = generation_of_the_row(n)
        mo_es = sum(var) / n
        dispersion_es = sum([(x - mo_es) ** 2 for x in var]) / (n - 1)
        interval_list.append((n - 1)*dispersion_es / chi2.isf((1 + (1 - a)) / 2, n - 1) -
                             (n - 1) * dispersion_es / chi2.isf((1 - (1 - a)) / 2, n - 1))
        dispersion_es_with_mo = sum([(x - mo) ** 2 for x in var]) / (n - 1)
        interval_list_with_mo.append(n*dispersion_es_with_mo / chi2.isf((1 + (1 - a)) / 2, n) -
                     n * dispersion_es_with_mo / chi2.isf((1 - (1 - a)) / 2, n))
    plt.plot(n_list, interval_list, label='Без известного матожидания')
    plt.plot(n_list, interval_list_with_mo, label='С известным МО')
    plt.title(
        "Зависимость величины доверительного интервала дисперсии\n от объёма выборки с доверительным значением " + str(1 - a))
    plt.legend()
    plt.show()


if __name__ == "__main__":
    n = 20
    # n = int(input("Enter n.\n"))
    variation = generation_of_the_row(n)
    mo_estimation = sum(variation) / n
    print("Точечная оценка матожидания: ", "%.4f" % mo_estimation)
    print("Матожидание:", "%.4f" % mo)
    dispersion_estimation = sum([(x - mo_estimation)**2 for x in variation]) / (n - 1)
    # dispersion_estimation = sum([(x**2 - mo_estimation**2) for x in variation]) / (n - 1)
    print("Точечная оценка дисперсии: ", "%.4f" % dispersion_estimation)
    print("Дисперсия - ",  "%.4f" % dispersion)
    s_estimated = dispersion_estimation ** 0.5
    confidence_intervals_for_expected_values(n, mo_estimation, s_estimated)
    confidence_interval_dependency(n, s_estimated)
    confidence_interval_volume_dependency()
    print()
    confidence_interval_for_all_dispersions(n, s_estimated, variation)
    confidence_interval_dependency_for_dispersion(n, variation, dispersion_estimation)
    confidence_interval_volume_dependency_for_dispersion()





