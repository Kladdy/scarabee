from enum import Enum

class Symmetry(Enum):
    """
    Defines the symmetry to use when simulating a PWR fuel assembly.
    """

    Full = 1
    """
    The fuel assembly has no symmetry.
    """

    Half = 2
    """
    The fuel assembly has symmetry about either the x or y axis (but not both).
    """

    Quarter = 3
    """
    The fuel assembly has symmetry about the x and y axis.
    """

