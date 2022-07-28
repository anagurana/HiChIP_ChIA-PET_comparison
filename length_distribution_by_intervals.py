import matplotlib.pyplot as plt

from length_freq_plots import frequency
from length_freq_plots import length
import math

global width
width = 7000


def freq_by_intervals(length, frequency, n_interval):
    # width = 7000
    k = 0
    coordinates = []

    interval = [width*n_interval, width*(n_interval+1)]

    for x in length:
        if interval[0] <= x <= interval[1]:
            coordinates.append(k)
        k = k + 1
    if not coordinates:
        print("There is no enought data")
    else:
        desired_freq = frequency[min(coordinates):max(coordinates)]
        return desired_freq, interval


def generate_histogram(n_interval):
    freq, interval = freq_by_intervals(length, frequency, n_interval)
    # plt.hist(freq, 1)
    plt.hist(freq, bins=range(min(freq), max(freq) + 1, 1), log=True)
    plt.title(f"Distribution of #PETs for {interval} bp length")
    plt.xlabel("Number of PETs")
    plt.ylabel("log(Frequency)")
    plt.show()


qnt_interval = math.ceil(max(length)/width)

generate_histogram(50)








