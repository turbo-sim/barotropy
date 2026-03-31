# Barotropy

Enhance your two-phase turbomachinery CFD simulations with **`barotropy`**, a Python package designed to create barotropic fluid property models. The barotropic approximation assumes that fluid properties (e.g., density, viscosity, speed of sound) depend solely on pressure, which is a very accurate assumption for most turbomachinery flows. This simplification reduces computational costs while improving convergence reliability in simulations involving complex fluid property variations, such as supercritical CO₂ compressors or two-phase flows in nozzles and turbines.

📦 **PyPI package**: [https://pypi.org/project/barotropy/](https://pypi.org/project/barotropy/)

📚 **Documentation**: [https://turbo-sim.github.io/barotropy/](https://turbo-sim.github.io/barotropy/) *(under construction)*  

🎓 **Tutorials**: [https://turbo-sim.github.io/barotropy/source/tutorials.html](https://turbo-sim.github.io/barotropy/source/tutorials.html)  
*Step-by-step examples for using the barotropic model in CFD solvers like Ansys Fluent and CFX.*


## Key Features

- **Simplified fluid modeling:** Generate barotropic models to lower computational cost and enhance solver robustness compared to real-gas fluid property tables.  
- **Suited for homogeneous two-phase flows:** Simulate two-phase flows with one or two-components using the Homogeneous Equilibrium Model (HEM) or the Delayed Equilibrium Model (DEM) assumptions..  
- **Seamless integration with commercial CFD solvers:** Easily export barotropic models as simple expressions ready to be copy-pasted into ANSYS Fluent or ANSYS CFX.  


<p align="center">
  <img src="docs/_static/two-phase_turbine.jpg"  width="85%" />
</p>

<p align="center">
  <img src="docs/_static/sCO2_compressor.jpg"  width="85%" />
</p>



## 🚀 User installation (via PyPI)

To begin using `barotropy`, install it via `pip`:

```bash
pip install barotropy
```

After installation, verify that everything is set up correctly by running the following command in your terminal:

```bash
python -c "import barotropy; barotropy.print_package_info()"
```

For detailed information and examples, visit the [documentation page](https://turbo-sim.github.io/barotropy/).


## 🛠️ Developer installation (from source)

To install `barotropy` from source, clone the repository and install it using [Poetry](https://python-poetry.org/):
```bash
git clone https://github.com/turbo-sim/barotropy.git
cd barotropy
poetry install
```

After installation, verify that everything is set up correctly by running:
```bash
poetry run python -c "import barotropy; barotropy.print_package_info()"
```


## License
The code in this repository is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.


## Contact Information

The code in this repository was developed by the [Sustainable Thermal Power group](https://thermalpower.dtu.dk/) at [DTU Construct](https://construct.dtu.dk/). Drop us an email at [roagr@dtu.dk](mailto:roagr@dtu.dk) if you have questions about the code or have a bug to report!

