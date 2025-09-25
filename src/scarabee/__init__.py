from ._scarabee import *
from . import _scarabee


# Check if we are in a Jupyter notebook or not
def _is_notebook() -> bool:
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False  # Probably standard Python interpreter


if _is_notebook():
    _scarabee._use_python_sink()

__author__ = _scarabee.__author__
__copyright__ = _scarabee.__copyright__
__license__ = _scarabee.__license__
__maintainer__ = _scarabee.__maintainer__
__email__ = _scarabee.__email__
__version__ = _scarabee.__version__
