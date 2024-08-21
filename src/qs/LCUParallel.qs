namespace LCUParallel {

    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;

    newtype CnotLayerSpec = (
    // specifies the stabilizers for the non-trivial resource state (z state / Omega_2 in the paper)
                        xstabs: Pauli[][],
                        zstabs: Pauli[][]
                        );


    newtype CliffordSpec = (
    // input from python preprocessing of required Clifford transformations
                            cnotL: CnotLayerSpec,
                            cnotM: CnotLayerSpec,
                            phase2Mask: Bool[],
                            cnotN: CnotLayerSpec,
                            phase1Mask: Bool[],
                            hadamardMask: Bool[],
                            // also include specs for inverse of the cnot layers
                            cnotLinv: CnotLayerSpec,
                            cnotMinv: CnotLayerSpec,
                            cnotNinv: CnotLayerSpec
                            );


    operation PrepareCSSFromStabs (resource_CSS: Qubit [], xstabs : Pauli [][], zstabs : Pauli [][] ) : Unit {
        // Prepares a CSS state from its stabilisers through measurement and correction
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
        // log depth fanout operation via recursion, can be converted to constant depth using ancilla or surface code compilation tricks
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

    operation CnotLayerPauliCorrection(qs: Qubit[], x_results: Result[], z_results: Result[], a_results: Result[]) : Unit {
        // Applies Pauli corrections after the constant depth CNOT layer based on the measurement results
        let n = Length(qs);
        // Working out the correction requires tracking the signs in the tableau, and correcting based on the measurement results
        // from the syndrome extraction, https://arxiv.org/abs/quant-ph/0406196, https://arxiv.org/abs/1805.12082, https://arxiv.org/abs/quant-ph/0408190
        // should be helpful for this, single Pauli layer not helpful for resource estimation so has been omitted for now
    }

    operation CnotLayerFromCnotLayerSpec(qs: Qubit[], cnotSpec: CnotLayerSpec): Unit {
        // Applies a CNOT circuit in constant depth using measurement and feedforward
        // CNOT layer specified by the stabilizers specifying the non-trivial resource state
        let n = Length(qs);
        use resource_ancillas = Qubit[2*n];
        use a_ancillas = Qubit[n];

        // cnot layer
        // prepare and use x resource state
        // very simple and can be done just by preparing bell pairs
        for i in 0 .. n - 1 {
            H(resource_ancillas[i]);
            CNOT(resource_ancillas[i], resource_ancillas[i+n]);
        }
        let x_results = HalfSteaneSyndromeMeasurement(qs, a_ancillas, resource_ancillas, true);
        ResetAll(resource_ancillas);


        // prepare and use z resource state
        PrepareCSSFromStabs(resource_ancillas, cnotSpec::xstabs, cnotSpec::zstabs);
        let z_results = HalfSteaneSyndromeMeasurement(qs, a_ancillas, resource_ancillas, false);

        let a_results = MeasureEachX(a_ancillas);

        CnotLayerPauliCorrection(qs, x_results, z_results, a_results);

        ResetAll(resource_ancillas);
        ResetAll(a_ancillas);
    }


    operation ApplyOnMask(qs: Qubit[], mask: Bool[], op: (Qubit => Unit)) : Unit {
        for i in 0 .. Length(qs) - 1 {
            if mask[i] {
                op(qs[i]);
            }
        }
    }

    operation CliffordConstantDepth(qs : Qubit [], clifford: CliffordSpec) : Unit is Adj {
        // Applies a Clifford circuit in constant depth using the canonical form
        // given in Proctor and Young https://arxiv.org/abs/2310.10882, applying CNOT layers in constant depth
        // using the scheme from https://arxiv.org/abs/1805.12082
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

        adjoint ... {
            ApplyOnMask(qs, clifford::hadamardMask, H);
            ApplyOnMask(qs, clifford::phase1Mask, Adjoint S);
            CnotLayerFromCnotLayerSpec(qs, clifford::cnotNinv);
            ApplyToEach(Adjoint S, qs);
            ApplyToEach(H, qs);
            ApplyOnMask(qs, clifford::phase2Mask, Adjoint S);
            CnotLayerFromCnotLayerSpec(qs, clifford::cnotMinv);
            ApplyToEach(Adjoint S, qs);
            CnotLayerFromCnotLayerSpec(qs, clifford::cnotLinv);
        }

    }


    operation ApplyLCUPartition(qs: Qubit[], prepRegisters: Qubit[][], nTerms: Int, clifford: CliffordSpec) : Unit is Adj {
        within {
            CliffordConstantDepth(qs, clifford);
        } apply {
            // apply hamiltonian terms (they have been transformed to PauliZs by the clifford circuit)
            for i in 0 .. nTerms - 1 {
                ApplyControlledOnInt(i, ApplyP(PauliZ, _), prepRegisters[i], qs[i]);
            }
        }
    }


    function SplitRegister(qs: Qubit[], k: Int) : Qubit[][] {
        // takes a register and splits it into k registers of equal size
        let n = Length(qs);
        let size = n / k;
        Fact(n % k == 0, "register must be divisible by k");
        mutable result: Qubit[][] = [];
        for i in 0 .. k-1 {
            set result += [qs[i*size .. (i+1)*size-1]];
        }
        return result;
    }

    @EntryPoint()
    operation testResource(): Unit {
        let n = 100;
        let nTerms = 1;
        use qs = Qubit[n];
        use ancillas = Qubit[nTerms * Ceiling(Lg(IntAsDouble(n)))];
        let prepRegisters = SplitRegister(ancillas, nTerms);
        for i in 0 .. nTerms - 1 {
            ApplyControlledOnInt(i, ApplyP(PauliZ, _), prepRegisters[i], qs[i]);
        }
    }
}

