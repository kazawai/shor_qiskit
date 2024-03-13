"""
RSA and Shor's algorithm.
Based on past Qiskit implementation : https://github.com/ttlion/ShorAlgQiskit/blob/master/Shor_Normal_QFT.py
"""

import traceback
# Time the execution of the algorithm
from time import time

from numpy import gcd, pi, zeros
from numpy.random import randint, seed
# For shor's algorithm in qiskit version 1.0.1
from qiskit import (ClassicalRegister, QuantumCircuit, QuantumRegister,
                    transpile)
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
# public_key, private_key = rsa(19391, 11939)
public_key, private_key = rsa(17, 11)


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


def qft(
    circuit: QuantumCircuit, qreg: QuantumRegister, n: int, with_swaps: bool = True
):
    """n-qubit QFT on q in circ."""
    for j in range(n):
        circuit.h(qreg[n - j - 1])
        for k in range(n - j - 1):
            circuit.cp(pi / float(2 ** (n - j - 1 - k)), qreg[k], qreg[n - j - 1])
        if with_swaps:
            circuit.barrier()

    if with_swaps:
        for j in range(n // 2):
            circuit.swap(qreg[j], qreg[n - j - 1])
    return circuit


def qft_dagger(
    circuit: QuantumCircuit, qreg: QuantumRegister, n: int, with_swaps: bool = True
):
    """n-qubit QFTdagger on q in circ."""
    if with_swaps:
        for j in range(n // 2):
            circuit.swap(qreg[j], qreg[n - j - 1])

    for j in range(n):
        for k in range(j + 1, n):
            circuit.cp(-pi / float(2 ** (k - j)), qreg[k], qreg[j])
        circuit.h(qreg[j])
        if with_swaps:
            circuit.barrier()
    return circuit


def get_angles(a, N):
    """Compute the angles for the QPE"""
    s = bin(int(a))[2:].zfill(N)
    angles = zeros([N])
    for i in range(0, N):
        for j in range(i, N):
            if s[j] == "1":
                angles[N - i - 1] += pow(2, -(j - i))
        angles[N - i - 1] *= pi
    return angles


def cc_phase_gate(
    circuit: QuantumCircuit,
    control_qubits: list[QuantumRegister],
    target_qubit: int,
    angle: float,
):
    """Create a doubly controlled phase gate"""
    circuit.cp(angle / 2, control_qubits[0], target_qubit)
    circuit.cx(control_qubits[1], control_qubits[0])
    circuit.cp(-angle / 2, control_qubits[0], target_qubit)
    circuit.cx(control_qubits[1], control_qubits[0])
    circuit.cp(angle / 2, control_qubits[1], target_qubit)


def phi_adder(
    circuit: QuantumCircuit,
    qreg: QuantumRegister,
    a: int,
    N: int,
    inverse: bool = False,
):
    """Create the circuit that performs the addition by a in the Fourier space (mod N)"""
    angles = get_angles(a, N)
    for i in range(N):
        if len(angles) == 0:
            break
        if not inverse:
            circuit.p(angles[i], qreg[i])
        else:
            circuit.p(-angles[i], qreg[i])


def controlled_phi_adder(
    circuit: QuantumCircuit,
    qreg: QuantumRegister,
    a: int,
    N: int,
    control_qubits: list[QuantumRegister] | QuantumRegister,
    inverse: bool = False,
):
    """Controlled version of the phi_adder"""
    angles = get_angles(a, N)
    for i in range(N):
        if len(angles) == 0:
            break
        if not inverse:
            if isinstance(control_qubits, QuantumRegister) or len(control_qubits) == 1:
                # Here, an error occurs when the control qubit and the target qubit are the same...
                try:
                    if control_qubits[0] != qreg[i]:
                        circuit.cp(
                            angles[i],
                            (
                                control_qubits
                                if isinstance(control_qubits, QuantumRegister)
                                else control_qubits[0]
                            ),
                            qreg[i],
                        )
                except Exception as e:
                    print(traceback.format_exc())
                    print(f"Control qubits: {control_qubits}")
                    print(f"Target qubit: {qreg[i]}")
                    print(f"Angle: {angles[i]}")
                    raise e

            else:
                cc_phase_gate(circuit, control_qubits, qreg[i], angles[i])
        else:
            if isinstance(control_qubits, QuantumRegister) or len(control_qubits) == 1:
                if control_qubits[0] != qreg[i]:
                    circuit.cp(
                        -angles[i],
                        (
                            control_qubits
                            if isinstance(control_qubits, QuantumRegister)
                            else control_qubits[0]
                        ),
                        qreg[i],
                    )
            else:
                cc_phase_gate(circuit, control_qubits, qreg[i], -angles[i])


def controlled_phi_adder_mod_N(
    circuit: QuantumCircuit,
    qreg: QuantumRegister,
    a: int,
    N: int,
    n: int,
    control_qubits: list[QuantumRegister] | QuantumRegister,
    ancilla_qubits: QuantumRegister,
):
    """Controlled version of the phi_adder"""
    controlled_phi_adder(circuit, qreg, a, n, control_qubits, False)
    phi_adder(circuit, qreg, N, n, True)
    qft_dagger(circuit, qreg, n, False)
    circuit.cx(qreg[n - 1], ancilla_qubits[0])
    qft(circuit, qreg, n, False)
    controlled_phi_adder(circuit, qreg, N, n, control_qubits, False)

    controlled_phi_adder(circuit, qreg, a, n, control_qubits, True)
    qft_dagger(circuit, qreg, n, False)
    circuit.x(qreg[n - 1])
    circuit.cx(qreg[n - 1], ancilla_qubits[0])
    circuit.x(qreg[n - 1])
    qft(circuit, qreg, n, False)
    controlled_phi_adder(circuit, qreg, a, n, control_qubits, False)


def controlled_phi_adder_mod_N_inv(
    circuit: QuantumCircuit,
    qreg: QuantumRegister,
    a: int,
    N: int,
    n: int,
    control_qubits: list[QuantumRegister] | QuantumRegister,
    ancilla_qubits: QuantumRegister,
):
    """Controlled version of the inverse phi_adder"""
    controlled_phi_adder(circuit, qreg, a, n, control_qubits, True)
    qft_dagger(circuit, qreg, n, False)
    circuit.x(qreg[n - 1])
    circuit.cx(qreg[n - 1], ancilla_qubits[0])
    circuit.x(qreg[n - 1])
    qft(circuit, qreg, n, False)
    controlled_phi_adder(circuit, qreg, a, n, control_qubits, False)

    controlled_phi_adder(circuit, qreg, N, n, control_qubits, False)
    qft_dagger(circuit, qreg, n, False)
    circuit.cx(qreg[n - 1], ancilla_qubits[0])
    qft(circuit, qreg, n, False)
    phi_adder(circuit, qreg, N, n, False)
    controlled_phi_adder(circuit, qreg, a, n, control_qubits, True)


def controller_mult_mod_N(
    circuit: QuantumCircuit,
    qreg: QuantumRegister,
    a: int,
    N: int,
    n: int,
    control_reg: list[QuantumRegister] | QuantumRegister,
    ancilla_reg: QuantumRegister,
):
    """Controlled version of the multiplication by a mod N"""
    qft(circuit, qreg, n, False)
    for i in range(n):
        controlled_phi_adder_mod_N(
            circuit, qreg, (2**i) * a % N, N, n, control_reg, ancilla_reg
        )
    qft_dagger(circuit, qreg, n, False)

    for i in range(n):
        if control_reg[0] != qreg[i]:
            circuit.cswap(control_reg[0], qreg[i], ancilla_reg[i])

    a_inv = pow(a, -1, N)
    qft(circuit, qreg, n + 1, False)

    for i in range(n):
        controlled_phi_adder_mod_N_inv(
            circuit,
            qreg,
            pow(2, i) * a_inv % N,
            N,
            n + 1,
            control_reg,
            ancilla_reg,
        )
    qft_dagger(circuit, qreg, n, False)


def check_power(N):
    """Check if N is a perfect power"""
    for i in range(2, N):
        if N % i == 0:
            if N % i**2 == 0:
                print(f"{N} is a perfect power")
                return i
            else:
                print(f"{N} is not a perfect power")
                break


def get_coprime(N):
    """Get a random coprime of N"""
    a = randint(2, N)
    while gcd(a, N) != 1:
        a = randint(2, N)
    return a


def qpe_period_finding(N, a, n) -> float | str:
    """Quantum Phase Estimation for period finding"""
    q_up = QuantumRegister(2 * n, "q_up")
    q_down = QuantumRegister(n, "q_down")
    a_q = QuantumRegister(n + 2, "a_q")
    c = ClassicalRegister(2 * n, "c")
    qpe = QuantumCircuit(q_up, q_down, a_q, c)
    try:
        # Initialize qubits
        qft(qpe, q_up, 2 * n, False)
        qpe.x(q_down[0])

        # Apply controlled multiplication gates
        for i in range(n):
            controller_mult_mod_N(qpe, q_up, a, N, n, q_down, a_q)

        # Apply inverse QFT
        qft_dagger(qpe, q_up, n, False)

        # Measure the qubits
        qpe.measure(q_up, c)

        # Transpile the circuit
        transpiled_qpe = transpile(qpe, AerSimulator())

        # Simulate the QuantumCircuit
        result = AerSimulator().run(transpiled_qpe).result()
        counts = result.get_counts(qpe)
        print(counts)

        # Extract the phase
        phase = max(counts, key=counts.get)
        phase = int(phase, 2) / 2**n

        print(qpe.draw())

        return phase
    except Exception as e:
        print(traceback.format_exc())
        # print(qpe.draw())
        return "Failure"


def shor(N):
    """Shor's algorithm"""
    # Check if N is a perfect power
    r = check_power(N)
    if r:
        return r, N // r

    # Get a random coprime of N
    a = get_coprime(N)

    # Check if a is a non-trivial square root of 1
    if pow(a, N // 2, N) == (N - 1):
        return gcd(a + 1, N), gcd(a - 1, N)

    # Find the period
    n = int(2 * (N**0.5)).bit_length()
    phase = qpe_period_finding(N, a, n)

    print(f"The phase is {phase}")

    # Check if the period is even
    if isinstance(phase, float) and phase % 1 == 0:
        return int((a ** (phase / 2) - 1) % N), int((a ** (phase / 2) + 1) % N)
    else:
        return "Failure"


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
