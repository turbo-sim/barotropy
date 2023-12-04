# Highlight exception messages
# https://stackoverflow.com/questions/25109105/how-to-colorize-the-output-of-python-errors-in-the-gnome-terminal/52797444#52797444
try:
    import IPython.core.ultratb
except ImportError:
    # No IPython. Use default exception printing.
    pass
else:
    import sys
    sys.excepthook = IPython.core.ultratb.ColorTB()


from glob import glob
from os.path import basename, dirname, isfile

# Get a list of all Python files in the current directory (excluding __init__.py)
modules = glob(dirname(__file__) + "/*.py")
__all__ = [basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]

# Import all modules specified in __all__
for module in __all__:
    exec(f"from .{module} import *")


# Import local functions

# from .fluid_properties import *
# from .fluent_automation import *
# from .fluent_plot_nozzle_data import *
# from .plot_options import *
# from .isentropic_nozzle import *
# from .PySolverView import *
# from .utilities import *
