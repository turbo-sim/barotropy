(developer_guide)=

# Developer guide

Thank you for considering contributing to this project! Here are some guidelines to help you get started:

## Installation from source

This installation guide is intended for developers who wish to contribute to or modify the `barotropy` source code. It assumes that the developer is using a Linux distribution or Windows with Git Bash terminal to have access to Git and Linux-like commands.

1. **Fork the repository**

   Go to the [project’s GitHub page](https://github.com/turbo-sim/barotropy) and click the **"Fork"** button in the upper right corner to create your own copy of the repository.

2. **Clone the forked repository**

   In your terminal, run:

   ```bash
   git clone https://github.com/<your-username>/barotropy.git
   cd barotropy
   ```

3. **Create a development environment**

   You can choose between two options:

   **Option A: Use Poetry directly**

   If you already have Poetry installed on your system, run:

   ```bash
   poetry install
   ```

   Poetry will automatically create a virtual environment and install the required dependencies.

   **Option B: Use Conda to manage the environment**

   If you prefer Conda for environment management, run:

   ```bash
   conda create --name barotropy_env python=3.13 pip poetry
   conda activate barotropy_env
   poetry install
   ```

   This sets up a Conda environment and installs Poetry inside it to manage dependencies.

4. **Verify the installation**

   Run the following to confirm the package is working:

   ```bash
   poetry run python -c "import barotropy; barotropy.print_package_info()"
   ```

   If the installation was successful, you will see the `barotropy` banner and package details printed in the terminal.

5. **Install additional packages**

   To add a runtime dependency:

   ```bash
   poetry add <package-name>
   ```

   To add a development-only dependency (e.g., for testing or docs):

   ```bash
   poetry add --dev <package-name>
   ```

   This will update both `pyproject.toml` and `poetry.lock` accordingly.

## Documentation

Install the development dependencies from the repository root:

```bash
poetry install --with dev
```

Start the local preview, which opens a browser and rebuilds when sources change:

```bash
poetry run python docs/build_docs.py
```

Build once (also used by CI):

```bash
poetry run python docs/build_docs.py --no-autobuild
```

The site is written to `docs/_build/html/index.html`. Builds fail on Sphinx
warnings so broken references and API imports are caught before deployment.
Use `--port 8080` to change the preview port, `--no-open-browser` to leave
the browser closed, or `--clean-only` to remove the HTML build. Paths are resolved relative to the script, so it also works from other
working directories when launched with the project's Python environment.

### Editing pages

Write pages in `docs/source/` as Markdown using MyST. Add new pages to a
`{toctree}` in `index.md` or the relevant section page. The `.rst` files in
`docs/source/api/` are generated from Python docstrings on every build; edit
the docstrings rather than these generated files. Private working notes stay
in `docs/dev_notes/`, outside the published source tree.

Sphinx settings, project paths, API exclusions, theme options, and Zotero
settings all live in `docs/conf.py`. The single `docs/build_docs.py` script
handles API generation, building, previewing, cleaning, and bibliography export.

### Citations and bibliography

Cite an entry with `{cite}` followed by its key in backticks, for example
{cite}`romei_computational_2021`. The bibliography page updates during both full
and incremental builds. Ordinary builds use the checked-in BibTeX file and do
not contact Zotero.

Refresh the bibliography from the configured Zotero group:

```bash
poetry run python docs/build_docs.py --build-bib --no-autobuild
```

Use `--build-bib-only` to export without building. Public groups need no API
key; for authenticated access, set the `ZOTERO_API_KEY` environment variable.
