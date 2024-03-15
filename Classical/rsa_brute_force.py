from time import time

from numpy import gcd, sqrt, ceil

N = 19391 * 11939
# N = 15


def classical_brute_force(N):
    for a in range(2, int(ceil(sqrt(N)))):
        if gcd(a, N) != 1:
            return a, N // a


if __name__ == "__main__":
    print("Algorithme classique de brute force")
    start = time()
    factors = classical_brute_force(N)
    print(N.bit_length())
    end = time()
    print(
        f"Les facteurs premiers de {N} sont {factors} : trouve en {end - start} secondes"
    )
    print("-" * 50)
