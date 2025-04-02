import os
import shutil
import subprocess
import argparse

def delete_build(build_dir="_build"):
    """
    Deletes previously generated documentation folders.

    Parameters:
    - build_dir (str): Directory where the built docs are stored.
    """
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)
        print(f"Deleted build directory: {build_dir}")

def run_sphinx_apidoc(output_dir, src_dir, exclude=None, force=False):
    """
    Run the sphinx-apidoc command to generate API documentation.

    Parameters:
    - source_dir (str): Directory where the source files reside.
    - output_dir (str): Directory to which the API documentation should be output.
    - module_path (str): Path to the module that should be documented.
    - extensions (list of str, optional): List of extensions to pass to the sphinx-apidoc command. 
                                          Default is ['-e'].

    Raises:
    - RuntimeError: If the sphinx-apidoc command fails.
    - FileNotFoundError: If the sphinx-apidoc command is not found.
    """
    cmd = ["sphinx-apidoc", "-o", output_dir, src_dir]

    # Exclude certain packages/modules
    if isinstance(exclude, list):
        for exclude_item in exclude:
            cmd.append(f"{src_dir}/{exclude_item}")

    cmd.append("-e") # put documentation for each module on its own page
    cmd.append('--no-toc') # Do not create a table of contents file
    cmd.append("-M") # Module first

    if force:
        cmd.append("-f") # Force overwriting on any existing generated files

    subprocess.check_call(cmd)
    print("Sphinx apidoc completed successfully.")


def run_sphinx_build(docs_dir=".", build_dir="_build", builder="html"):
    """
    Run sphinx-build to generate the documentation.

    Parameters:
    -----------
    - source_dir (str): Directory containing the .rst and other source files.
    - build_dir (str): Directory where the output will be written.
    - builder (str): The output format (e.g., "html", "latex").

    """
    cmd = ["sphinx-build", "-b", builder, docs_dir, build_dir]
    subprocess.check_call(cmd)
    print(f"Sphinx build ({builder} format) completed successfully.")

def run_sphinx_autobuild(docs_dir=".", build_dir="_build", port=8000):
    """
    Run sphinx-autobuild to serve and auto-rebuild docs on changes.

    Parameters:
    - docs_dir (str): Directory with documentation source files.
    - build_dir (str): Output directory.
    - port (int): Port to serve HTML.
    """
    cmd = [
        "sphinx-autobuild",
        docs_dir,
        build_dir,
        "--open-browser",
        "--port", str(port),
    ]
    subprocess.check_call(cmd)



def run_script(script_path):
    """Run the nomenclature generation script."""
    cmd = ["python", script_path]
    subprocess.check_call(cmd)



EXCLUDE_MODULES = [
    "../barotropy/mcerp",
    "../barotropy/matlab",
    "../barotropy/uncertainty_utilities.py",
    "../barotropy/kaist_postprocessing.py",
]

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Build or preview Sphinx documentation.")
    parser.add_argument("--no-autobuild", action="store_true", help="Disable live preview (sphinx-autobuild) and use sphinx-build instead.")
    args = parser.parse_args()

    delete_build()
    # run_script("build_nomenclature.py")
    # run_script("build_bibliography.py")
    # run_script("build_configuration.py")
    if args.no_autobuild:
        run_sphinx_apidoc(output_dir="source/api/", src_dir="../barotropy", exclude=EXCLUDE_MODULES)
        run_sphinx_build(docs_dir=".", build_dir="_build")
    else:
        run_sphinx_apidoc(output_dir="source/api/", src_dir="../barotropy", exclude=EXCLUDE_MODULES)
        run_sphinx_autobuild(docs_dir=".", build_dir="_build")
