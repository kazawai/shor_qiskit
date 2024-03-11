from time import time

from matplotlib import pyplot as plt

from rsa_brute_force import classical_brute_force
from rsa_shor import shor

# Run the algorithm for different values of N
primes = [631, 641, 2297, 2309, 4001, 4003, 10007, 10009, 11939, 19391]
times = []
values = [p * q for p in primes for q in primes if p != q]
for value in values:
    start = time()
    shor(value)
    end = time()
    times.append(end - start)
c_times = []
for value in values:
    start = time()
    classical_brute_force(value)
    end = time()
    c_times.append(end - start)

# Plot the results in a graph with same scale on the y-axis
plt.plot(values, times, label="Quantum")
plt.plot(values, c_times, label="Classical")
plt.xlabel("Value of N")
plt.ylabel("Time (s)")
plt.legend()
plt.title("Shor's algorithm vs brute force")

plt.show()
