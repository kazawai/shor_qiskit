# Time the execution of the algorithm
from time import time

from numpy import gcd, pi
from numpy.random import seed
# For shor's algorithm in qiskit version 1.0.1
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator


def rsa(P, Q):
    """Very simple RSA key generation function."""
    N = P * Q
    phi = (P - 1) * (Q - 1)
    E = 3
    while gcd(E, phi) != 1:
        E += 2
    D = 0
    while D * E % phi != 1:
        D += 1
    return (E, N), (D, N)


# Define message
message = "HELLO"
# Generate the public and private keys
public_key, private_key = rsa(19391, 11939)


def rsa_encrypt(message, public_key):
    """Very simple RSA encryption function."""
    E, N = public_key
    return [pow(ord(char), E, N) for char in message]


def rsa_decrypt(message, private_key):
    """Very simple RSA decryption function."""
    D, N = private_key
    return "".join([chr(pow(char, D, N)) for char in message])


# Print the encrypted message
ciphertext = rsa_encrypt(message, public_key)
print(f"Encrypted text : {[hex(c) for c in ciphertext]}")

# Print the decrypted message
decrypted = rsa_decrypt(ciphertext, private_key)
print(f"Decrypted text : {decrypted}")

print("-" * 50)

seed(1)

N = public_key[1]

print("La cle publique est", public_key)
print("-" * 50)

####################
# Shor's algorithm #
####################


def c_amod15(a, power):
    """Controlled multiplication by a mod 15"""
    if a not in [2, 7, 8, 11, 13]:
        raise ValueError("'a' must be 2,7,8,11, or 13")
    U = QuantumCircuit(4)
    for _ in range(power):
        if a in [2, 13]:
            U.swap(0, 1)
            U.swap(1, 2)
            U.swap(2, 3)
        if a in [7, 8]:
            U.swap(2, 3)
            U.swap(1, 2)
            U.swap(0, 1)
        if a == 11:
            U.swap(1, 3)
            U.swap(0, 2)
        if a in [7, 11, 13]:
            for q in range(4):
                U.x(q)
    U = U.to_gate()
    U.name = "%i^%i mod 15" % (a, power)
    c_U = U.control()
    return c_U


def qft_dagger(n):
    """n-qubit QFTdagger the first n qubits"""
    qc = QuantumCircuit(n)
    for qubit in range(n // 2):
        qc.swap(qubit, n - 1 - qubit)
    for j in range(n):
        for m in range(j):
            qc.cp(-pi / float(2 ** (j - m)), m, j)
        qc.h(j)
    qc.name = "QFT†"
    return qc


def qpe_period_finding():
    """Estimate the period of the function a^x mod 15"""
    a = 7
    n_count = 8
    qc = QuantumCircuit(4 + n_count, n_count)

    # Apply Hadamard gates
    for q in range(n_count):
        qc.h(q)  # Initialize counting qubits in state |+>

    qc.x(3 + n_count)  # And auxiliary register in state |1>
    # Apply controlled-U operations
    for q in range(n_count):  # Do controlled-U operations
        qc.append(c_amod15(a, 2**q), [q] + [i + 4 for i in range(n_count)])

    # Apply the inverse-QFT
    qc.append(qft_dagger(n_count), range(n_count))
    qc.measure(range(n_count), range(n_count))

    aer_sim = AerSimulator()

    t_qc = transpile(qc, aer_sim)
    qobj = aer_sim.run(t_qc).result()

    counts = qobj.get_counts()
    return counts


def denominator(r: float):
    return 1 / r - int(1 / r)


# Define the Shor's algorithm function
def shor(N):
    """Shor's Algorithm"""
    if N % 2 == 0:
        return 2, N // 2

    # Check if N is a perfect power
    for a in range(2, N):
        if N % a == 0:
            if pow(a, N - 1, N) == 1:
                continue
            else:
                return a, N // a
    # Call the quantum part of the algorithm
    counts = qpe_period_finding()

    # Process the results
    for measured_value in counts:
        print(f"Measured {int(measured_value, 2)}")
        phase = int(measured_value, 2) / (2**8)
        frac = phase - int(phase)
        print(f"Phase {phase} = {frac} (mod 1)")
        if frac < 0.5:
            r = frac
        else:
            r = frac - 1
        print(f"Closest fraction {r}")

        # Check if the period is even
        if denominator(r) % 2 == 0:
            print(f"Period {r}")
            x = a ** (denominator(r) // 2)
            print(f"x = {a}^{denominator(r)//2} = {x} (mod {N})")
            if (x + 1) % N != 0:
                print(f"Non-trivial factor {gcd(x + 1, N)}")
                return gcd(x + 1, N), gcd(x - 1, N)
    return None


if __name__ == "__main__":
    # Time the execution of the algorithm
    print("Algorithme quantique de Shor")
    start = time()
    factors = shor(N)
    end = time()
    print(
        f"Les facteurs premiers de {N} sont {factors} : trouve en {end - start} secondes"
    )
    print("-" * 50)
