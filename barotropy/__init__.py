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
from .fluent_automation import *
from .fluid_properties import *
from .barotropic_model import *

# Import submodules
from .math import *
from .plot_options import *
from .sCO2_utilities import *
from .utilities import *

# Package info
__version__ = "0.0.2"
PACKAGE_NAME = "barotropy"
URL_GITHUB = "https://github.com/turbo-sim/barotropy"
URL_DOCS = "https://turbo-sim.github.io/barotropy/"
URL_DTU = "https://thermalpower.dtu.dk/"
BREAKLINE = 80 * "-"

def print_banner():
    """Prints a banner."""
    banner = """
                ____                   __                        
               / __ )____ __________  / /__________  ____  __  __
              / __  / __ `/ ___/ __ \/ __/ ___/ __ \/ __ \/ / / /
             / /_/ / /_/ / /  / /_/ / /_/ /  / /_/ / /_/ / /_/ / 
            /_____/\__,_/_/   \____/\__/_/   \____/ .___/\__, /  
                                                 /_/    /____/   
"""
#     banner = """
#             ____  ___    ____  ____  __________  ____  ______  __
#            / __ )/   |  / __ \/ __ \/_  __/ __ \/ __ \/ __ \ \/ /
#           / __  / /| | / /_/ / / / / / / / /_/ / / / / /_/ /\  / 
#          / /_/ / ___ |/ _, _/ /_/ / / / / _, _/ /_/ / ____/ / /  
#         /_____/_/  |_/_/ |_|\____/ /_/ /_/ |_|\____/_/     /_/   
                                                                
# """
    print(BREAKLINE)
    print(banner)
    print(BREAKLINE)
    # https://manytools.org/hacker-tools/ascii-banner/
    # Style: Sland


def print_package_info():
    """Prints package information with predefined values."""

    info = f""" Version:       {__version__}
 Repository:    {URL_GITHUB}
 Documentation: {URL_DOCS}"""
    print_banner()
    print(info)
    print(BREAKLINE)




