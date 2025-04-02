.. _developer_guide:

Developer guide
=======================

Thank you for considering contributing to this project! Here are some guidelines to help you get started:


Installation from source
----------------------------

This installation guide is intended for developers who wish to contribute to or modify the ``barotropy`` source code. It assumes that the developer is using a Linux distribution or Windows with Git Bash terminal to have access to Git and Linux-like commands.

1. **Fork the repository**

   Go to the `project’s GitHub page <https://github.com/turbo-sim/barotropy>`_ and click the **"Fork"** button in the upper right corner to create your own copy of the repository.


2. **Clone the forked repository**

   In your terminal, run:

   .. code-block:: bash

      git clone https://github.com/<your-username>/barotropy.git
      cd barotropy


3. **Create a development environment**

   You can choose between two options:

   **Option A: Use Poetry directly**

   If you already have Poetry installed on your system, run:

   .. code-block:: bash

      poetry install

   Poetry will automatically create a virtual environment and install the required dependencies.

   **Option B: Use Conda to manage the environment**

   If you prefer Conda for environment management, run:

   .. code-block:: bash

      conda create --name barotropy_env python=3.13 pip poetry
      conda activate barotropy_env
      poetry install

   This sets up a Conda environment and installs Poetry inside it to manage dependencies.

4. **Verify the installation**

   Run the following to confirm the package is working:

   .. code-block:: bash

      poetry run python -c "import barotropy; barotropy.print_package_info()"

   If the installation was successful, you will see the ``barotropy`` banner and package details printed in the terminal.


5. **Install additional packages**

   To add a runtime dependency:

   .. code-block:: bash

      poetry add <package-name>

   To add a development-only dependency (e.g., for testing or docs):

   .. code-block:: bash

      poetry add --dev <package-name>

   This will update both `pyproject.toml` and `poetry.lock` accordingly.


.. 
   Pull request guidelines
   -------------------------

   Please follow these steps to submit a pull request.

   1. **Create a branch in your forked repository**:

      - Open your terminal in the projects root.
      - Create branch:

      .. code-block:: bash

         git checkout -b <feature-name>

   2. **Make your changes**:

      - Implement your feature or bugfix.


   3. **Commit your changes**:

      .. code-block:: bash 

         git commit -m "Description of changes"

   4. **Push to your fork**: 

      .. code-block:: bash

         git push origin feature-name

   5. **Open a pull request**: 

      - Go to your fork on GitHub and click the "New pull request" button.


.. 
   Reporting issue
   ----------------

   If you find a bug or have a feature request, please open an issue in the Github project page and follow the provided templates.

   CI/CD Pipeline
   --------------

   barotropy uses GitHub Actions to automate its Continuous Integration and Continuous Deployment (CI/CD) processes.

   Automated Testing
   ^^^^^^^^^^^^^^^^^

   The ``ci.yml`` action is triggered whenever a commit is pushed to the repository. This action runs the test suite on both Windows and Linux environments, ensuring the code's compatibility and correctness across different platforms.

   Package Publishing
   ^^^^^^^^^^^^^^^^^^

   barotropy utilizes the ``bumpversion`` package to manage versioning and release control. To increment the version number, use the following command:

   .. code-block:: bash

      bumpversion patch  # or minor, major

   After bumping the version, push the changes to the remote repository along with tags to signify the new version:

   .. code-block:: bash

      git push origin --tags

   If the tests pass successfully, the package is automatically published to the Python Package Index (PyPI), making it readily available for users to install and use.

   Documentation Deployment
   ^^^^^^^^^^^^^^^^^^^^^^^^

   barotropy automates the deployment of documentation using the ``deploy_docs`` action. This action builds the Sphinx documentation of the project and publishes the HTML files to GitHub Pages each time that a new commit is pushed to the remote repository. By automating this process, barotropy ensures that the project's documentation remains up-to-date and easily accessible to users and contributors.