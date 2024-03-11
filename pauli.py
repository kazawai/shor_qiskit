from qiskit import (ClassicalRegister, QuantumCircuit, QuantumRegister,
                    transpile)
from qiskit_aer import AerSimulator

# Define the quantum circuit
qc = QuantumCircuit(QuantumRegister(1), ClassicalRegister(1))

# Append the X operator
qc.x(0)

# Measure the qubit
qc.measure(0, 0)

# Display the circuit
print(qc)

# Transpile the circuit for the AerSimulator
aer_sim = AerSimulator()

aer_sim_transpile = transpile(qc, aer_sim)

# Simulate the transpiled circuit 10 times
result = aer_sim.run(aer_sim_transpile, shots=10).result()

# Print the result
print(result.get_counts())
