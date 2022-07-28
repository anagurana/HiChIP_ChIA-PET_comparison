import matplotlib.pyplot as plt

from length_freq_plots import frequency as frequencies
from length_freq_plots import length as lengths
import math

loop_len = 10000


def freq_by_deslength(length, frequency, loop_length):
    k = 0
    coordinates = []
    if loop_length >= 1000000:
        for x in length:
            if x <= loop_length:
                coordinates.append(k)
            k = k + 1
        else:
            for x in length:
                if (loop_length - 1000) <= x <= (loop_length + 1000):
                    coordinates.append(k)
                k = k + 1
    desired_freq = frequency[min(coordinates):max(coordinates)]
    return desired_freq


def generate_histogram(loop_length):
    freq = freq_by_deslength(lengths, frequencies, loop_length)
    plt.hist(freq, bins=range(min(freq), max(freq) + 1, 1))
    plt.title(f"Distribution of frequency for {loop_length} bp length")
    plt.xlabel('frequency (#PETs)')
    plt.ylabel('quantity')
    plt.show()


# qnt_interval = math.ceil(max(lengths)/width)

generate_histogram(loop_len)
