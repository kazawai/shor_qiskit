from numpy import pi, random
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit_aer import AerSimulator

backend = AerSimulator()


def bell_state():
    q = QuantumRegister(2)
    c = ClassicalRegister(2)
    circuit = QuantumCircuit(q, c)

    circuit.h(q[0])
    circuit.cx(q[0], q[1])
    return circuit


def algorithm(value1, value2, circuit):
    theta1, theta2 = 0, pi / 8
    if value1 == 1:
        theta1 = pi / 4
    if value2 == 1:
        theta2 = -theta2

    circuit.ry(theta1, 0)
    circuit.ry(theta2, 1)

    return circuit


def quantum_strategy(x, y):
    circuit = bell_state()
    circuit = algorithm(x, y, circuit)
    circuit.measure([0, 1], [0, 1])

    result = backend.run(circuit).result()
    counts = result.get_counts(circuit)

    counts = dict(counts)

    return int(list(counts.keys())[0][0]), int(list(counts.keys())[0][1])


def bell_inequality_game(nb_tries):
    wins = 0
    for _ in range(nb_tries):
        x = random.randint(2)
        y = random.randint(2)
        a, b = quantum_strategy(x, y)
        wins += x * y == (a and b)
    print(f"Quantum strategy wins {wins}/{nb_tries} times")


if __name__ == "__main__":
    bell_inequality_game(1000)
