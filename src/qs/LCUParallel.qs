namespace LCUParallel {

    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;

    operation PrepareCSSFromStabs (resource_CSS: Qubit [], xstabs : Pauli [][], zstabs : Pauli [][] ) : Unit {
        Fact(Length(xstabs) == Length(zstabs), "xstabs and zstabs must have the same length");
        Fact(Length(xstabs) == Length(resource_CSS), "must have same number of stabilizers and qubits");

        // reset qubits to state 0, Z stabilizers all correct
        // ResetAll(resource_CSS);

        // measure all X stabilizers
        mutable results: Result [] = [];
        for xstab in xstabs {
            set results += [Measure(xstab, resource_CSS)];
        }

        // Fix the state by applying the appropriate Z stabilizers
        for index in 0 .. Length(results)-1 {
            if results[index] == One {
                ApplyPauli(zstabs[index], resource_CSS);
            }
        }

    }


    operation CNOTTransversal(reg1: Qubit[], reg2: Qubit[]) : Unit is Adj + Ctl {
        Fact(Length(reg1) == Length(reg2), "registers must have the same length");
        for i in 0 .. Length(reg1) - 1 {
            CNOT(reg1[i], reg2[i]);
        }
    }


    operation FanoutRegister(origin_register: Qubit[], target_registers: Qubit[][]) : Unit is Adj + Ctl {
        // log depth fanout operation via recursion, can be converted to linear depth using ancilla or surface code compilation tricks
        // TEST THIS
        let n_reg = Length(target_registers);
        if (n_reg == 1) {
            CNOTTransversal(origin_register, target_registers[0]);
        }

        if (n_reg >= 2) {
            CNOTTransversal(origin_register, target_registers[n_reg/2]);
            FanoutRegister(origin_register, target_registers[0..n_reg/2]);
            FanoutRegister(target_registers[n_reg/2], target_registers[n_reg/2+1...]);
        }
    }

    operation MeasureEachX(qs: Qubit[]) : Result[] {
        ApplyToEach(H, qs);
        return MeasureEachZ(qs);
    }


    operation HalfSteaneSyndromeMeasurement(qs: Qubit[], ancilla : Qubit [], resource_CSS : Qubit [], x_type: Bool) : Result [] {
        // Takes a register qs on n qubits and a resource state on 2n qubits and performs steane syndrome extraction, using another n ancilla
        // Can be done in parallel using 2 2n resource states

        mutable results: Result [] = [];
        let n = Length(qs);
        Fact(Length(resource_CSS) == 2*n, "resource state must be 2 times the size of the target register");
        Fact(Length(ancilla) == n, "ancilla must be the same size as the target register");
        

        // x measurements
        if (x_type) {
            CNOTTransversal(resource_CSS[...n], ancilla);
            CNOTTransversal(resource_CSS[n...], ancilla);
            set results += MeasureEachX(resource_CSS[0..n]);
        } else { // z measurements
            CNOTTransversal(ancilla, resource_CSS[...n]);
            CNOTTransversal(qs, resource_CSS[n...]);
            set results += MeasureEachZ(resource_CSS[n...]);
        }
    
        return results;
    }

    operation CnotLayerFromCnotLayerSpec(qs: Qubit[], cnotSpec: CnotLayerSpec): Unit {
        let n = Length(qs);
        use resource_ancillas = Qubit[2*n];
        use a_ancillas = Qubit[n];

        // cnot layer
        // prepare and use x resource state
        PrepareCSSFromStabs(resource_ancillas, cnotSpec::xstabs1, cnotSpec::zstabs1);
        let x_results = HalfSteaneSyndromeMeasurement(qs, a_ancillas, resource_ancillas, true);
        // prepare and use z resource state
        PrepareCSSFromStabs(resource_ancillas, cnotSpec::xstabs2, cnotSpec::zstabs2);
        let z_results = HalfSteaneSyndromeMeasurement(qs, a_ancillas, resource_ancillas, false);

        let a_results = MeasureEachX(a_ancillas);

        // do pauli corrections (SHOULD BE POSSIBLE TO FIGURE OUT)


        ResetAll(resource_ancillas);
        ResetAll(a_ancillas);
    }

    newtype CnotLayerSpec = (xstabs1: Pauli[][], 
                            zstabs1: Pauli[][], 
                            xstabs2: Pauli[][],
                            zstabs2: Pauli[][]
                            );

    newtype CliffordSpec = (cnotL: CnotLayerSpec,
                            cnotM: CnotLayerSpec,
                            phase2Mask: Bool[],
                            cnotN: CnotLayerSpec,
                            phase1Mask: Bool[],
                            hadamardMask: Bool[]
                            );


    operation ApplyOnMask(qs: Qubit[], mask: Bool[], op: (Qubit => Unit)) : Unit {
        for i in 0 .. Length(qs) - 1 {
            if mask[i] {
                op(qs[i]);
            }
        }
    }

    operation CliffordConstantDepth(qs : Qubit [], clifford: CliffordSpec) : Unit is Adj {
        // docstring
        body ... {      
            CnotLayerFromCnotLayerSpec(qs, clifford::cnotL);
            ApplyToEach(S, qs);
            CnotLayerFromCnotLayerSpec(qs, clifford::cnotM);
            ApplyOnMask(qs, clifford::phase2Mask, S);
            ApplyToEach(H, qs);
            ApplyToEach(S, qs);
            CnotLayerFromCnotLayerSpec(qs, clifford::cnotN);
            ApplyOnMask(qs, clifford::phase1Mask, S);
            ApplyOnMask(qs, clifford::hadamardMask, H);       
        }

        // think about how to do adjoint version of this, just need to correct resource states
        adjoint ... {
            // calculate the adjoint clifford spec... etC
        }

    }

    operation ApplyLCUPartition(qs: Qubit[], prepRegisters: Qubit[][], hamilTerms, clifford: CliffordSpec) : Unit is Adj {
        body ... {
            // Apply clifford transformation on main register
            CliffordConstantDepth(qs, clifford);
            // apply hamiltonian terms
            for i in 0 .. Length(hamilTerms) - 1 {
                let term = hamilTerms[i];
                let reg = prepRegisters[i];
                ApplyControlledOnInt(i, termFunc, reg, qs);
            }
            // apply clifford adjoint
            Adjoint CliffordConstantDepth(qs, clifford);
        }
    }

}


