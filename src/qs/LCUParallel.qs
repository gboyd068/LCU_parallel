namespace LCUParallel {

    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;

    operation PrepareCSSFromStabs (qs: Qubit [], xstabs : Pauli [][], zstabs : Pauli [][] ) : Unit {
        Fact(Length(xstabs) == Length(zstabs), "xstabs and zstabs must have the same length");
        Fact(Length(xstabs) == Length(qs), "must have same number of stabilizers and qubits");

        // reset qubits to state 0, Z stabilizers all correct
        // ResetAll(qs);

        // measure all X stabilizers
        mutable results: Result [] = [];
        for xstab in xstabs {
            set results += [Measure(xstab, qs)];
        }

        // Fix the state by applying the appropriate Z stabilizers
        for index in 0 .. Length(results)-1 {
            if results[index] == One {
                ApplyPauli(zstabs[index], qs);
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
        // currently linear depth fanout, use log depth! or constant depth that zhu showed us for surface code
        let n = Length(target_registers);
        for i in 0 .. n - 1 {
            CNOTTransversal(origin_register, target_registers[i]);
        }
    }


    operation ApplyResourceState (qs: Qubit[], resource_CSS : Qubit [] ) : Unit is Adj + Ctl {
        // Do steane(?) syndrome extraction


        // Do the correction
    }

    operation CliffordConstantDepth(qs : Qubit []) : Unit is Adj + Ctl {
        // think about how to do adjoint and controlled versions of this
        body ... {
            // Do the magic
            let n = Length(qs);
            use ancillas = Qubit[4*n];

            // cnot layer


            // ancillas go out of scope here
            ResetAll(ancillas);
        }

        controlled (ctls, ...) { 

        }

        adjoint ... {
            
        }
    }
}