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

# Import local functions

from .fluid_properties import *
from .fluent_automation import *
from .fluent_plot_nozzle_data import *
from .plot_options import *
from .isentropic_nozzle import *
from .PySolverView import *
from .utilities import *