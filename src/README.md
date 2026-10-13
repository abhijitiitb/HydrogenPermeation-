# NEB Configuration Generation for NiH and PdH

## Overview

This repository contains Fortran 90 programs and LAMMPS input files for generating initial and final atomic configurations and performing Nudged Elastic Band (NEB) calculations for hydrogen migration in nickel (NiH) and palladium (PdH) systems.

The workflow generates hydrogen occupation configurations, constructs the corresponding atomic structures, relaxes the initial and final configurations using LAMMPS, and prepares the endpoints for NEB calculations.

## Repository Structure

```text
NEB-Million-Configurations/
├── Ni/
│   ├── NEB_OCCUP1.txt
│   ├── NiAlH_jea.eam.alloy
│   ├── excute1_4000Nih.f90
│   ├── in_NIH_initial.eam
│   ├── in_NIH_initial1.eam
│   ├── input_NEB
│   ├── occupation_in_1NNN1.txt
│   ├── print_occupation_4000_files_final.f90
│   └── resultP1.sh
├── Pd/
│   ├── NEB_OCCUP1.txt
│   ├── PdAgH_HybridPd3Ag.eam.alloy
│   ├── excute1_4000PdH.f90
│   ├── in_PdH_initial.eam
│   ├── in_PdH_initial1.eam
│   ├── input_NEB
│   ├── occupation_in_1NNN1.txt
│   ├── print_occupation_4000_files_final.f90
│   └── resultP1.sh
├── .gitignore
└── README.md
```

## Directory Descriptions

### Ni — Hydrogen Migration in Nickel

The `Ni/` directory contains the Fortran source programs, occupation files, EAM/alloy potential, and LAMMPS input files required for the NiH workflow.

### Pd — Hydrogen Migration in Palladium

The `Pd/` directory contains the corresponding Fortran source programs, occupation files, EAM/alloy potential, and LAMMPS input files for the PdH workflow.

## File Descriptions

| File | Purpose |
|---|---|
| `NEB_OCCUP1.txt` | Contains occupation configurations to process. |
| `occupation_in_1NNN1.txt` | Stores an occupation vector used to construct a particular configuration. |
| `excute1_4000Nih.f90` | Fortran source program for the NiH configuration-generation workflow. |
| `excute1_4000PdH.f90` | Fortran source program for the PdH configuration-generation workflow. |
| `print_occupation_4000_files_final.f90` | Fortran source for generating initial and final hydrogen occupation vectors. |
| `NiAlH_jea.eam.alloy` | Interatomic potential file used by the NiH workflow. |
| `PdAgH_HybridPd3Ag.eam.alloy` | Interatomic potential file used by the PdH workflow. |
| `in_NIH_initial.eam` | LAMMPS input for relaxing the initial NiH configuration. |
| `in_NIH_initial1.eam` | LAMMPS input for relaxing the final NiH configuration. |
| `in_PdH_initial.eam` | LAMMPS input for relaxing the initial PdH configuration. |
| `in_PdH_initial1.eam` | LAMMPS input for relaxing the final PdH configuration. |
| `input_NEB` | LAMMPS input for the NEB calculation. |
| `resultP1.sh` | Shell script for automating the configuration-processing workflow. |

## Workflow

The general workflow for each material is:

1. **Generate occupation configurations:** Read occupation patterns from `NEB_OCCUP1.txt`.
2. **Construct occupation vectors:** Use the Fortran occupation-generation program to create initial and final hydrogen occupation vectors.
3. **Generate atomic structures:** Use the material-specific Fortran program to construct the initial and final configurations.
4. **Relax the structures:** Run the corresponding LAMMPS input files to relax the initial and final configurations.
5. **Prepare NEB endpoints:** Extract the required endpoint structures and migration information.
6. **Perform the NEB calculation:** Use `input_NEB` to calculate the energy profile along the hydrogen migration pathway.
7. **Collect results:** Use `resultP1.sh` to automate processing and collect the relevant output data.

## Software Requirements

- Fortran compiler, such as GNU Fortran (`gfortran`).
- LAMMPS compiled with the required EAM/alloy potential support and NEB functionality.
- MPI implementation if running LAMMPS in parallel.
- Bash, `awk`, `sed`, and standard Unix utilities.

## Running the Workflow

Run the workflow separately within the `Ni/` and `Pd/` directories.

The exact compilation and execution commands depend on the compiler, LAMMPS executable, MPI configuration, and filenames specified in the source programs and shell scripts.

Before running a complete production workflow:

1. Verify that the potential file is compatible with the LAMMPS `pair_style` and `pair_coeff` commands.
2. Compile the Fortran programs.
3. Confirm that the input and output filenames match across the Fortran, shell, and LAMMPS files.
4. Test a small number of configurations before generating a large dataset.
5. Verify the relaxed endpoint structures and NEB energy profiles.

## Reproducibility

The Fortran source files, occupation inputs, potential files, and LAMMPS input scripts are provided to document the configuration-generation and NEB workflows.

Simulation results depend on the interatomic potential, simulation parameters, initial configurations, and software environment.

## License and Potential Files

The licensing terms for the source code and interatomic potential files should be checked before redistribution or reuse. Include the appropriate attribution and license information where required.
