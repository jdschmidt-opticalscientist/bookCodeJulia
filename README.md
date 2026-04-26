# bookCodeJulia: Numerical Simulation of Optical Wave Propagation

This repository provides the official Julia implementation of the simulation framework from the textbook **"Numerical Simulation of Optical Wave Propagation with Examples in MATLAB"**.

It transitions the core wave optics algorithms—including Fresnel and Fraunhofer propagation, phase screen generation, and atmospheric turbulence modeling—into an accessible Julia environment.

---

## 📚 Reference & Citation

If you use this code in your research or professional work, please cite the original textbook:

> Jason D. Schmidt, _Numerical Simulation of Optical Wave Propagation with Examples in MATLAB_, SPIE, Bellingham, WA (2010).

- [SPIE Digital Library](https://www.spiedigitallibrary.org/ebooks/PM/Numerical-Simulation-of-Optical-Wave-Propagation-with-Examples-in-MATLAB/eISBN-9780819483270/10.1117/3.866274)
- [Amazon Reference](https://www.amazon.com/Numerical-Simulation-Propagation-Examples-Monograph/dp/0819483265/)


---

## 📝 A Note on the Julia Implementation

**Direct Correspondence:** The primary goal of this repository is to maintain a one-to-one correspondence with the MATLAB examples provided in the textbook.

To assist readers in connecting the book material with the Julia implementation, the code is intentionally written to mirror the MATLAB structure, including folder structure, variable names, and comments. As a result, it is not optimized for maximum execution speed, nor does it strictly follow "Julia" idioms that might obscure the underlying physics for those following along with the book.

**Enhancement:** Unlike the original code provided on the disc with the book, this Julia version includes integrated plotting routines to visualize key outputs immediately.

---

## 📂 Project Structure
The core package and its environment are contained within the `OpticalWavePropSim` subdirectory:

```text
bookCodeJulia/
├── examples/           # Textbook example scripts
├── OpticalWavePropSim/ # Julia package environment
│   ├── src/            # Core propagation source code
│   ├── Project.toml    # Dependencies
│   └── Manifest.toml   # Locked environment state
```

---

## 🚀 Installation & Quick Start

### 1. Requirements

This project requires Julia 1.9 or later.

### 2. Setup

Clone the repository and instantiate the project environment from the root directory. This ensures all required packages (FFTW, JuMP, Ipopt, Plots) are installed in a reproducible state:

```bash
git clone [https://github.com/jdschmidt-opticalscientist/bookCodeJulia.git](https://github.com/jdschmidt-opticalscientist/bookCodeJulia.git)
cd bookCodeJulia/OpticalWavePropSim
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 3. Run an Example

To run any of the simulation examples, execute them using the project environment flag (--project) pointed at the OpticalWavePropSim directory:

```julia
# Run from the repository root:
julia --project=OpticalWavePropSim examples/Chapter9/example_pt_source_atmos_setup.jl
```

---

##⚖️ License

The Julia source code in this repository is licensed under the BSD 3-Clause License. See the LICENSE file for full details.

Note: This license applies specifically to this Julia translation. The original MATLAB code and textbook content remain under the copyright of SPIE.

## Development Quality

![Julia CI](https://github.com/jdschmidt-opticalscientist/bookCodeJulia/actions/workflows/main.yml/badge.svg?branch=main)

We use GitHub Actions to ensure code quality through automated linting and unit testing.
