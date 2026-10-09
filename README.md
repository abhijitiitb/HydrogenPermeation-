# Activation Energy Dataset

The activation energies are evaluated using CI-NEB method at discrete hydrogen concentrations (x) (fraction of hydrogen occupied octahedral sites) in NiH_x and PdH_x systems. The mapping of activation energy barriers of metal hydride system, hydrogen concentration with their corresponding file is listed in the table below.

| Metal Name | Hydrogen Concentration  |             Filename             |
|:----------:|:-----------------------:|:--------------------------------:|
|   NiH_x    |          0.001          | NiH_activation_energy_x0.001.txt |
|   NiH_x    |          0.25           |  NiH_activation_energy_x025.txt  |
|   NiH_x    |          0.50           | NiH_activation_energy_x0.50.txt  |
|   NiH_x    |          0.75           | NiH_activation_energy_x0.75.txt  |
|   NiH_x    |          0.99           | NiH_activation_energy_x0.99.txt  |
|   PdH_x    |          0.001          | PdH_activation_energy_x0.001.txt |
|   PdH_x    |          0.25           | PdH_activation_energy_x0.25.txt  |
|   PdH_x    |          0.50           | PdH_activation_energy_x0.50.txt  |
|   PdH_x    |          0.99           | PdH_activation_energy_x0.99.txt  |

Each file consists of 262,144 (=2^18) energy barriers for all 1NN configurations of 18 sites, where each local environment is represented by an 18-bit binary configuration corresponding to the occupancy of the 18 1NN sites:

- Site 1 represents the Most Significant Bit (MSB).
- Site 18 represents the Least Significant Bit (LSB).
- 1 indicates that the corresponding site is occupied by hydrogen.
- 0 indicates that the corresponding site is unoccupied.
- Conversion of the binary configuration to decimal yields the 0-based row index, with the first row corresponding to index 0.
- The unit of the activation energy barrier is eV.

## Example

If we want to find the activation energy barrier for NiH_x at hydrogen concentration of 0.25 with out of 18 sites if only site ID \[2,4,5\] are occupied then will get the index of the row corresponding to this local environment configuration in binary as follows,

| Site ID  | 1   | 2    | 3  | 4   | 5   | 6  | 7  | ... | 18 |
|----------|-----|------|----|-----|-----|----|----|-----|----|
| Occupied | No  | Yes  | No | Yes | Yes | No | No | ... | No |
| 0/1      | 0   | 1    | 0  | 1   | 1   | 0  | 0  | ... | 0  |

Binary Configuration = 010110000000000000 

Converting From Binary to Decimal, we get 90112

Row Number = 90112 (starting with first row as index 0) 

or

Row Number = 90113 (starting with first row as index 1)

Thus checking the corresponding row number of text file NiH_activation_energy_x025.txt, we get

Activation energy barrier = 0.5292291 eV