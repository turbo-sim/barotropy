(developer_guide)=

# Developer guide

Thank you for considering contributing to this project! Here are some guidelines to help you get started:

## Installation from source

With Git and Poetry installed, clone the repository and install its dependencies:

```bash
git clone https://github.com/turbo-sim/barotropy.git
cd barotropy
poetry install
```

Poetry manages the development environment and installs the package and its dependencies.
If you plan to submit changes through a fork, use your fork's URL in the clone command.

Verify the installation:

```bash
poetry run python -c "import barotropy; barotropy.print_package_info()"
```

### Adding dependencies

To add a runtime dependency:

```bash
poetry add <package-name>
```

To add a development dependency:

```bash
poetry add --group dev <package-name>
```

These commands update both `pyproject.toml` and `poetry.lock`.

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
the browser closed, or `--clean-only` to remove the HTML build. Paths are resolved relative to the script, so it also works from other working directories when launched with the project's Python environment.

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
