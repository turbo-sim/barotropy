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

# Import subpackages
from .fluent import *
from .props import *
from .barotropic_model import *

# Import submodules
from .plot_options import *
from .sCO2_utilities import *
from .utilities import *

# Set plot options
set_plot_options(minor_ticks=False)

