from pysolver_view import *

# Import subpackages
from .properties import *
from .fluent_automation import *

# Import submodules
from .barotropic_model import *
from .math import *
from .graphics import *
from .utilities import *

# Other code (to be cleaned)
from .sCO2_utilities import *
from .isentropic_nozzle import *

# Package info
__version__ = "0.1.1"
PACKAGE_NAME = "barotropy"
URL_GITHUB = "https://github.com/turbo-sim/barotropy"
URL_DOCS = "https://turbo-sim.github.io/barotropy/"
URL_DTU = "https://thermalpower.dtu.dk/"
BREAKLINE = 80 * "-"

def print_banner():
    """Prints a banner."""
    banner = r"""
                ____                   __                        
               / __ )____ __________  / /__________  ____  __  __
              / __  / __ `/ ___/ __ \/ __/ ___/ __ \/ __ \/ / / /
             / /_/ / /_/ / /  / /_/ / /_/ /  / /_/ / /_/ / /_/ / 
            /_____/\__,_/_/   \____/\__/_/   \____/ .___/\__, /  
                                                 /_/    /____/   
"""
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





