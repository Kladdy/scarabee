from .._scarabee import FormFactors
from ..reseau import Symmetry
import numpy as np


class CoreFormFactors:
    """
    Evlauates the power form factors for any position across the entire core of
    the reactor. It is used to accomplish pin-power reconstruction after full
    core simulations.

    Parameters
    ----------
    tile_width : float
        Width of a single tile. Generally the assembly width/pitch.
    num_tiles : int
        Number of tiles across the entire core at the widest point.
    form_factors : 3D array of FormFactors or 0.
        The array contains the form factors for all the assemblies in the core
        at all axial slices. It should be ordered in what you see is what you
        get (z,y,x), and the indexing for z and y is backwards (low z index
        means high z value). Entries for reflectors or zero flux regions should
        contain a 0.
    z_widths : 1D array of float
        Width of each axial slice from low to high z.
    symmetry : Symmetry
        Symmetry used for the simulation. Default value is Symmetry.Full.
    """

    def __init__(
        self,
        tile_width: float,
        num_tiles: int,
        form_factors: np.ndarray,
        z_widths: np.ndarray,
        symmetry: Symmetry = Symmetry.Full,
    ):
        if tile_width <= 0.0:
            raise ValueError("Tile width must be > 0.")
        if num_tiles <= 0.0:
            raise ValueError("Number of tiles must be > 0.")
        if form_factors.ndim != 3:
            return TypeError("Array of form factors must have 3 dimensions.")
        if z_widths.ndim != 1:
            return TypeError("Array of z widths must have 1 dimension.")
        if np.min(z_widths) <= 0.0:
            raise ValueError("All widths must be > 0.")
        if z_widths.size != form_factors.shape[0]:
            raise IndexError(
                "The number of z widths does not agree with the number of axial slices."
            )
        if symmetry == Symmetry.Half:
            raise ValueError("Only Full or Quarter symmetry are currently supported.")

        self.tile_width = tile_width
        self.num_tiles = num_tiles
        self.form_factors = form_factors
        self.symmetry = symmetry

        # Make sure that the shape agrees with the symmetry !
        expected_num_tiles = num_tiles
        if symmetry == Symmetry.Quarter:
            expected_num_tiles = num_tiles // 2 + (num_tiles % 2)
        if self.form_factors.shape[1] != expected_num_tiles:
            raise RuntimeError(
                "Shape of form factors along y does not agree with symmetry and number of tiles."
            )
        if self.form_factors.shape[2] != expected_num_tiles:
            raise RuntimeError(
                "Shape of form factors along x does not agree with symmetry and number of tiles."
            )

        # Create the widths arrays
        self.z_bounds = np.zeros(z_widths.size + 1)
        for k in range(1, self.z_bounds.size):
            self.z_bounds[k] = self.z_bounds[k - 1] + z_widths[k - 1]

        self.x_bounds = np.zeros(expected_num_tiles + 1)
        self.y_bounds = np.zeros(expected_num_tiles + 1)
        if self.symmetry == Symmetry.Quarter:
            self.x_bounds[1] = 0.5 * self.tile_width
            self.y_bounds[1] = 0.5 * self.tile_width
        else:
            self.x_bounds[1] = self.tile_width
            self.y_bounds[1] = self.tile_width
        for i in range(2, self.x_bounds.size):
            self.x_bounds[i] = self.x_bounds[i - 1] + self.tile_width
            self.y_bounds[i] = self.y_bounds[i - 1] + self.tile_width

        # Check scalars and floats >= 0.
        for k in range(self.form_factors.shape[0]):
            for j in range(self.form_factors.shape[1]):
                for i in range(self.form_factors.shape[2]):
                    if isinstance(self.form_factors[k, j, i], int):
                        self.form_factors[k, j, i] = float(self.form_factors[k, j, i])
                    if (
                        isinstance(self.form_factors[k, j, i], float)
                        and self.form_factors[k, j, i] < 0.0
                    ):
                        raise ValueError("All form factors must be >= 0.")

        # Make sure the widths agree with the symmetry and tile width
        for k in range(self.form_factors.shape[0]):
            for j in range(self.form_factors.shape[1]):
                for i in range(self.form_factors.shape[2]):
                    tile = self.form_factors[-(k + 1), -(j + 1), i]

                    if isinstance(tile, float):
                        continue

                    if (
                        i == 0
                        and self.symmetry == Symmetry.Quarter
                        and abs(tile.x_width - 0.5 * self.tile_width) > 1.0e-6
                    ):
                        raise ValueError(
                            "Width along x does not agree with symmetry and tile width."
                        )
                    elif i == 0 and abs(tile.x_width - self.tile_width) > 1.0e-6:
                        raise ValueError(
                            "Width along x does not agree with symmetry and tile width."
                        )

                    if (
                        j == -self.form_factors.shape[1]
                        and self.symmetry == Symmetry.Quarter
                        and abs(tile.y_width - 0.5 * self.tile_width) > 1.0e-6
                    ):
                        raise ValueError(
                            "Width along y does not agree with symmetry and tile width."
                        )
                    elif (
                        j == -self.form_factors.shape[1]
                        and abs(tile.y_width - self.tile_width) > 1.0e-6
                    ):
                        raise ValueError(
                            "Width along y does not agree with symmetry and tile width."
                        )

    def __call__(self, x, y, z):
        """
        Evaluates the form factors at the provided points.

        Parameters
        ----------
        x : float or np.ndarray of float
            x coordinate at which to evaluate the form factors
        y : float or np.ndarray of float
            y coordinate at which to evaluate the form factors
        z : float or np.ndarray of float
            z coordinate at which to evaluate the form factors

        Returns
        -------
        float or np.ndarray
            If x, y, and z are all floats, then a single float with a form
            factor is returned. If any of them are an array, an array of form
            factors is returned.
        """
        # Check if any argument was an array
        has_array = (
            isinstance(x, np.ndarray)
            or isinstance(y, np.ndarray)
            or isinstance(z, np.ndarray)
        )

        # Convert float to array if needed
        if has_array:
            if not isinstance(x, np.ndarray):
                x = np.array([float(x)])
            if not isinstance(y, np.ndarray):
                y = np.array([float(y)])
            if not isinstance(z, np.ndarray):
                z = np.array([float(z)])

        if not has_array:
            return self._ff_single_point(x, y, z)

        # Make sure all arrays are 1D
        if x.ndim != 1:
            raise ValueError("Array of x values must be 1D.")
        if y.ndim != 1:
            raise ValueError("Array of y values must be 1D.")
        if z.ndim != 1:
            raise ValueError("Array of z values must be 1D.")

        # Initialize output array
        out = np.zeros((x.size, y.size, z.size))

        # Evaluate the form factor at all points
        for i in range(x.size):
            for j in range(y.size):
                for k in range(z.size):
                    out[i, j, k] = self._ff_single_point(x[i], y[j], z[k])
        return out

    def _ff_single_point(self, x, y, z) -> float:
        # We first find the Z slice
        k = 0
        while z > self.z_bounds[k + 1]:
            k += 1
            if k == self.z_bounds.size - 1:
                raise ValueError(
                    f"The provided z value of {z} is above the maximum z value of {self.z_bounds[-1]}."
                )

        # Find the y index in true space
        j = 0
        while y > self.y_bounds[j + 1]:
            j += 1
            if j == self.y_bounds.size - 1:
                raise ValueError(
                    f"The provided y value of {y} is above the maximum y value of {self.y_bounds[-1]}."
                )

        # Find the x index in true space
        i = 0
        while x > self.x_bounds[i + 1]:
            i += 1
            if i == self.x_bounds.size - 1:
                raise ValueError(
                    f"The provided x value of {x} is above the maximum x value of {self.x_bounds[-1]}."
                )

        try:
            ff = self.form_factors[-(k + 1), -(j + 1), i]
            if isinstance(ff, FormFactors):
                return ff(x - self.x_bounds[i], y - self.y_bounds[j])
            else:
                return ff
        except:
            raise RuntimeError(
                f"Could not evaluate form factor at ({x}, {y}) translated to ({x - self.x_bounds[i]}, {y - self.y_bounds[j]})."
            )
