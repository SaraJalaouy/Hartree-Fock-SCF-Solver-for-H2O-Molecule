
# 💧 Hartree-Fock SCF Solver for H₂O Molecule

This project implements a self-consistent field (SCF) solver based on the **Hartree-Fock method** for the **water molecule (H₂O)** using **Fortran**. It computes the molecular orbitals, density matrix, Fock matrix, Hartree-Fock energy, and MP2 correlation correction from input integral data.

## 🧪 Features

- Computes overlap, one-electron, and two-electron integrals
- Builds and diagonalizes the Fock matrix
- Performs SCF iterations with convergence checking
- Outputs Hartree-Fock energy
- Transforms molecular orbital coefficients
- Calculates MP2 (second-order Møller-Plesset) energy correction
- Outputs in Molden format
- Computes atomic charges

## 📂 Project Structure
  ├── main.f90 # Main program file
  ├── hfmod.f90 # Module containing supporting subroutines and functions
  ├── wasser.out # Output file in Molden format 


  ## 📥 Input Files

The program expects precomputed integrals for the H₂O molecule:
- Overlap matrix
- One-electron integrals
- Two-electron integrals
- Nuclear repulsion energy

## 🚀 How to Run

1. Compile the program using a Fortran compiler:

```bash
gfortran -o hfh2o main.f90 hfmod.f90 matrixmulti.f90 -L libprak.a
````
2. Run the executable
  hfh2o (program_name)


📘 Output
- SCF iteration and convergence information
- Final Hartree-Fock energy
- MP2 energy correction
- Molden-compatible output file
- Atomic charges


