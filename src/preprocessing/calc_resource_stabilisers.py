import stim

# test required stim functionality

t = stim.Tableau(3)
H = stim.Tableau.from_named_gate("H")
cnot = stim.Tableau.from_named_gate("CNOT")
t.append(H, [0])
t.append(cnot, [0, 1])

t = stim.Tableau(3)
identity3 = stim.Tableau.from_conjugated_generators(
        xs=[
            stim.PauliString("X__"),
            stim.PauliString("_X_"),
            stim.PauliString("__X"),
        ],
        zs=[
            stim.PauliString("Z__"),
            stim.PauliString("_Z_"),
            stim.PauliString("__Z"),
        ],
    )

t = stim.Tableau.from_stabilizers(stabilizers=[
            stim.PauliString("XY_"),
            stim.PauliString("_YX"),
            stim.PauliString("ZZZ"),],
            allow_underconstrained=True)
# The conjugation of P by C is equal to C**-1 * P * C.

print(t.z_output(2)) 

