# Excipient Effects of Arginine on Hydrophobic Interactions and Biomolecular Models
Thank you for your interest in our project! This repository contains all necessary files and scripts to perform molecular simulations and subsequent analyses for our hydrophobic polymer / arginine, lysozyme / arginine, and porcine parvovirus / arginine systems. Below, you'll find an overview of the repository structure and how to get started.

## Citation
If you use this repository in your research, please cite the following publications:

> Jonathan W. P. Zajac, Praveen Muralikrishnan, Idris Tohidian, Xianci Zeng, Caryn L. Heldt, Sarah L. Perry, Sapna Sarupria. 16, 6780-6792, (2025) Chemical Science. DOI: 10.1039/d4sc08672d

> Jonathan W. P. Zajac, Praveen Muralikrishnan, Caryn L. Heldt, Sarah L. Perry, Sapna Sarupria. Advance Article, (2025) Molecular Systems Design & Engineering. DOI: 10.1039/D4ME00201F

## Repository Structure

### **Simulations** (`simulations/`)
This directory contains all files required to set up and run molecular dynamics (MD) simulations. It is organized into the following subdirectories:

- **`hydrophobic-polymer/`** - All scripts and input required to reproduce our hydrophobic polymer simulations.
- **`hewl/`** - Additional input and parameters required to reproduce lysozyme simulations.
- **`ppv/`** - Additional input and parameters required to reproduce porcine parvovirus simulations.
- **`mdp/`** – Stores MD parameter files (`.mdp`) defining simulation settings.
- **`mdrun/`** – Scripts and configurations for executing MD runs.
- **`topol/`** – Topology files (`.top`, `.itp`) describing molecular structures and interactions.
- **`box-solv/`** – Scripts and files for defining the simulation box and solvating the system.
- **`insert-molecules/`** – Scripts for inserting excipient molecules into the simulation box.
- **`gen-ions/`** – Configuration files for adding ions to neutralize the system.
- **`reus/`** - Scripts and tools used to run REUS simulations via PLUMED.
- **`ndx/`** - All index files generated to aid in energy group definitions or analysis.
- **`conf/`** - Initial configurations used for simulations.

### **Analysis** (`analysis/`)
This directory contains scripts and tools for post-simulation analysis, including structure validation, energy calculations, and trajectory processing.

## Getting Started
To use this repository:
1. Clone the repository:
   ```bash
   git clone git@github.com:SAMPEL-Group/Arginine-Hydrophobic-Polymer.git
   cd Arginine-Hydrophobic-Polymer
   ```
2. Ensure dependencies (e.g., GROMACS, Python libraries) are installed as needed.
3. Download the CHARMM22 force field from the Alex MacKerell website.
4. Download the TIP4P/2005 water model from the Carlos Vega website.
5. Merge force field files provided in `simulations/hydrophobic-polymer/topol` with the CHARMM force field directory.
        - Move `hydr-poly.rtp` to `charmm.ff/` directory.
        - Paste the contents of `hydr-poly-ffnonbonded.itp` to `charmm.ff/ffnonbonded.itp`.
        - Paste the contents of `hydr-poly-ffbonded.itp` to `charmm.ff/ffbonded.itp`.
        - Paste the contents of `hydr-poly-atomtypes.itp` to `charmm.ff/atomtypes.atb`.
        - Place `tip4p05.itp` in the `charmm.ff/` directory.

## Contributions
Contributions are welcome! Please submit a pull request or open an issue if you have suggestions or improvements.

## License
All written and graphical materials here are made available under a CC-BY 4.0 license, and all source code/software is made available under an MIT license. Both of these allow broad reuse with attribution.

## Contact
For questions or feedback, please contact Jonathan Zajac at zajac028@umn.edu.
