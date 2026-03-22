from .._scarabee import (
    DiffusionData,
    FormFactors,
)
from typing import Optional
from abc import ABC, abstractmethod
import copy


class CoreTile(ABC):
    """
    A CoreTile represents a unique geometry unit in a a nodal full core
    simulation, typically a single fuel assembly or reflector region.
    """

    def __init__(self):
        pass

    @property
    @abstractmethod
    def num_x_slots(self) -> int:
        """
        Number of unique nodes which should be used along the x direction.
        """
        return 0

    @property
    @abstractmethod
    def num_y_slots(self) -> int:
        """
        Number of unique nodes which should be used along the y direction.
        """
        return 0

    @property
    @abstractmethod
    def x_width(self) -> float:
        """
        Width of the tile along the x direction.
        """
        return 0.0

    @property
    @abstractmethod
    def y_width(self) -> float:
        """
        Width of the tile along the y direction.
        """
        return 0.0

    @abstractmethod
    def rotate_clockwise(self):
        """
        Rotates the ADFs and form factors, corresponding to a 90 degree
        rotation of the assembly in the clockwise direction.
        """
        return self

    @abstractmethod
    def rotate_counterclockwise(self):
        """
        Rotates the ADFs and form factors, corresponding to a 90 degree
        rotation of the assembly in the counter-clockwise direction.
        """
        return self

    @abstractmethod
    def reflect_across_x_axis(self):
        """
        Performs a reflection of the ADFs and form factors across the x axis.
        This means that the +y and -y ADFs are swapped.
        """
        return self

    @abstractmethod
    def reflect_across_y_axis(self):
        """
        Performs a reflection of the ADFs and form factors across the y axis.
        This means that the +x and -x ADFs are swapped.
        """
        return self


class SimpleTile(CoreTile):
    """
    Represents a single fuel assembly or reflector region, with diffusion cross
    sections, ADFs, CDFs, and form factors. Only a single set of cross sections
    are provided for the entire assembly. This is appropriate for symmetric
    fuel assemblies or reflector region.

    Parameters
    ----------
    diffusion_data : DiffusionData
        The diffusion cross sections, ADFs, and CDFs for the tile.
    form_factors : FormFactors
        Form factors for the tile.
    """

    def __init__(self, diffusion_data: DiffusionData, form_factors: FormFactors):
        self.diffusion_data = diffusion_data
        self.form_factors = form_factors

    @property
    def num_x_slots(self) -> int:
        return 1

    @property
    def num_y_slots(self) -> int:
        return 1

    @property
    def x_width(self) -> float:
        return self.form_factors.x_width

    @property
    def y_width(self) -> float:
        return self.form_factors.y_width

    def rotate_clockwise(self) -> CoreTile:
        self.diffusion_data.rotate_clockwise()
        self.form_factors.rotate_clockwise()
        return self

    def rotate_counterclockwise(self) -> CoreTile:
        self.diffusion_data.rotate_counterclockwise()
        self.form_factors.rotate_counterclockwise()
        return self

    def reflect_across_x_axis(self) -> CoreTile:
        self.diffusion_data.reflect_across_x_axis()
        self.form_factors.reflect_across_x_axis()
        return self

    def reflect_across_y_axis(self) -> CoreTile:
        self.diffusion_data.reflect_across_y_axis()
        self.form_factors.reflect_across_y_axis()
        return self


class QuadrantsTile(CoreTile):
    """
    Represents a single fuel assembly with diffusion cross sections, ADFs, CDFs,
    and form factors. Each quadrant of the assembly has a unique DiffusionData
    instance. This is appropriate for asymmetric fuel assemblies. A single
    FormFactors instance is provided, which can be generated from independent
    quadrant assembly calculations.

    Parameters
    ----------
    quad1 : DiffusionData or None
        The diffusion cross sections, ADFs, and CDFs for the I quadrant.
    quad2 : DiffusionData or None
        The diffusion cross sections, ADFs, and CDFs for the II quadrant.
    quad3 : DiffusionData or None
        The diffusion cross sections, ADFs, and CDFs for the III quadrant.
    quad4 : DiffusionData or None
        The diffusion cross sections, ADFs, and CDFs for the IV quadrant.
    form_factors : FormFactors
        Form factors for the tile.
    """

    def __init__(
        self,
        quad1: Optional[DiffusionData],
        quad2: Optional[DiffusionData],
        quad3: Optional[DiffusionData],
        quad4: Optional[DiffusionData],
        form_factors: FormFactors,
    ):
        non_none_count = len([q for q in [quad1, quad2, quad3, quad4] if q is not None])

        if non_none_count == 3:
            raise RuntimeError(
                "Cannot define a QuadrantsTile with 3 quadrants. Must provide 1, 2, or 4 quadrants."
            )

        if non_none_count == 2 and (
            (quad1 is None and quad3 is None) or (quad2 is None and quad4 is None)
        ):
            raise RuntimeError(
                "If providing only 2 quadrants to a QuadrantsTile, they cannot be diagonal to one another."
            )

        self.quad1 = quad1
        self.quad2 = quad2
        self.quad3 = quad3
        self.quad4 = quad4
        self.form_factors = form_factors

    @staticmethod
    def from_independent_quadrant(
        diffusion_data: DiffusionData, form_factors: FormFactors
    ) -> "QuadrantsTile":
        """
        Creates a QuadrantsTile from a single DiffusionData and FormFactors
        instance which were generated as independent quadrants, where ADFs and
        CDFs were generated along the symmetry axes of the assembly.

        Parameters
        ----------
        diffusion_data : DiffusionData
            Few-group cross sections, ADFs, and CDFs for the quarter assembly.
        form_factors : FormFactors
            Form factors for the quarter assembly.

        Returns
        -------
        QuadrantsTile
        """
        q1_dd = copy.deepcopy(diffusion_data)
        q1_ff = copy.deepcopy(form_factors)
        q2_dd = copy.deepcopy(q1_dd)
        q2_ff = copy.deepcopy(q1_ff)
        q2_dd.rotate_counterclockwise()
        q2_ff.rotate_counterclockwise()
        q3_dd = copy.deepcopy(q2_dd)
        q3_ff = copy.deepcopy(q2_ff)
        q3_dd.rotate_counterclockwise()
        q3_ff.rotate_counterclockwise()
        q4_dd = copy.deepcopy(q3_dd)
        q4_ff = copy.deepcopy(q3_ff)
        q4_dd.rotate_counterclockwise()
        q4_ff.rotate_counterclockwise()
        return QuadrantsTile(
            q1_dd, q2_dd, q3_dd, q4_dd, FormFactors(q1_ff, q2_ff, q3_ff, q4_ff)
        )

    @property
    def num_x_slots(self) -> int:
        if (self.quad1 is not None and self.quad2 is not None) or (
            self.quad3 is not None and self.quad4 is not None
        ):
            return 2
        else:
            return 1

    @property
    def num_y_slots(self) -> int:
        if (self.quad1 is not None and self.quad4 is not None) or (
            self.quad2 is not None and self.quad3 is not None
        ):
            return 2
        else:
            return 1

    @property
    def x_width(self) -> float:
        return self.form_factors.x_width

    @property
    def y_width(self) -> float:
        return self.form_factors.y_width

    def rotate_clockwise(self) -> CoreTile:
        tmp = [self.quad1, self.quad2, self.quad3, self.quad4]
        for dd in tmp:
            if dd is not None:
                dd.rotate_clockwise()
        self.quad1 = tmp[1]
        self.quad2 = tmp[2]
        self.quad3 = tmp[3]
        self.quad4 = tmp[0]
        self.form_factors.rotate_clockwise()
        return self

    def rotate_counterclockwise(self) -> CoreTile:
        tmp = [self.quad1, self.quad2, self.quad3, self.quad4]
        for dd in tmp:
            if dd is not None:
                dd.rotate_counterclockwise()
        self.quad1 = tmp[3]
        self.quad2 = tmp[0]
        self.quad3 = tmp[1]
        self.quad4 = tmp[2]
        self.form_factors.rotate_counterclockwise()
        return self

    def reflect_across_x_axis(self) -> CoreTile:
        tmp = [self.quad1, self.quad2, self.quad3, self.quad4]
        for dd in tmp:
            if dd is not None:
                dd.reflect_across_x_axis()
        self.quad1 = tmp[3]
        self.quad2 = tmp[2]
        self.quad3 = tmp[1]
        self.quad4 = tmp[0]
        self.form_factors.reflect_across_x_axis()
        return self

    def reflect_across_y_axis(self) -> CoreTile:
        tmp = [self.quad1, self.quad2, self.quad3, self.quad4]
        for dd in tmp:
            if dd is not None:
                dd.reflect_across_y_axis()
        self.quad1 = tmp[1]
        self.quad2 = tmp[0]
        self.quad3 = tmp[3]
        self.quad4 = tmp[2]
        self.form_factors.reflect_across_y_axis()
        return self
