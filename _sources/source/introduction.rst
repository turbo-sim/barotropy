
.. _installation:



Getting started
===============

User installation guide
------------------------

This guide shows how to install ``barotropy`` using ``pip``. The package can be installed in your system Python, but we recommend using a separate **virtual environment** to avoid conflicts with other packages.

You can choose between **Conda** or the built-in Python ``venv`` module depending on your preference.

1. **Basic installation with pip**

   If you already manage your environments and want a quick install, run:

   .. code-block:: bash

      pip install barotropy

   This will install the latest version of ``barotropy`` from PyPI.

2. **Recommended: install in an isolated environment**

   Choose one of the following options to create a clean environment and install ``barotropy``:

   **Option 1: Using Conda**

   If you use Conda or Miniconda to manage environments:

   .. code-block:: bash

      conda create --name barotropy_env python=3.13 -y
      conda activate barotropy_env
      pip install barotropy

   **Option 2: Using venv (lightweight)**

   If you prefer to avoid Conda and use Python’s built-in ``venv`` module:

   - On Windows (with Git Bash or PowerShell):

     .. code-block:: bash

        python -m venv env
        source env/Scripts/activate
        python -m pip install barotropy

   - On Linux or macOS:

     .. code-block:: bash

        python3 -m venv env
        source env/bin/activate
        python3 -m pip install barotropy

3. **Verify the installation**

   After installation, check that ``barotropy`` is correctly installed:

   .. code-block:: bash

      python -c "import barotropy; barotropy.print_package_info()"

   If successful, the console will display a banner and package information.

Congratulations! You're ready to start using ``barotropy`` to generate fluid property models.




Minimal working example
-----------------------

To be completed

