# Barotropy


This package contains functions for generating a barotropic model for fluid properties, a simplified mathematical representation that assumes fluid properties are only a function of pressure. This approach can be particularly useful in CFD modeling, as it offers a faster and more robust alternative to more complex fluid models. Barotropic models are particularly relevant for certain types of flows, such as two-phase flowsin turbomachinery. In such applications, the barotropic model is advantageous because it can approximate scenarios where the fluid follows a known thermodynamic process, such as an adiabatic process within a turbine or compressor.


## Getting Started

To begin using `barotropy`, install it via `pip`:

```bash
pip install barotropy
```

After installation, verify that everything is set up correctly by running the following command in your terminal:

```bash
python -c "import barotropy; barotropy.print_package_info()"
```

For detailed information and examples, visit the [documentation page](https://turbo-sim.github.io/barotropy/) and explore the [tutorials](https://turbo-sim.github.io/barotropy/source/tutorials.html) to get started with the code!



## License
The code in this repository is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.


## Contact Information

The code in this repository was developed by the [Sustainable Thermal Power group](https://thermalpower.dtu.dk/) at [DTU Construct](https://construct.dtu.dk/). Drop us an email at [roagr@dtu.dk](mailto:roagr@dtu.dk) if you have questions about the code or have a bug to report!

