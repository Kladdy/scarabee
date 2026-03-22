from .._scarabee import (
    Vector,
    Direction,
    MOCDriver,
    BoundaryCondition,
    CMFD,
    ADF,
    CDF,
)
from .nodal_flux import NodalFlux1D, NodalFlux2D
from .symmetry import Symmetry
import numpy as np
from typing import List, Tuple


def compute_adf_cdf_from_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    symmetry: Symmetry,
    independent_quadrant: bool,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the assembly and corner discontinuity factors in the few-group
    structure using the CMFD results. Currently, this method neglects the
    coupled x-y terms in the nodal construction when computing the CDFs.
    This simplifies the algorithm, but also ensures problems won't occur if
    trying to run single-pin cell "assembly" calculations.

    Parameters
    ----------
    moc : MOCDriver
        Solved MOC problem for which the equivalence factors will be computed.
    moc_to_nodal_cond_scheme : list of list of ints
        Scheme used to condense from the MOC problem to the resulting nodal
        diffusion parameters.
    symmetry : Symmetry
        Indicates the symmetry of the problem.
    independent_quadrant : bool
        If quarter symmetry is indicated (symmetry = Symmetry.Quarter), then
        unique ADFs and CDFs will still be generated for each side, as if full
        symmetry were being used.

    Returns
    -------
    ADF : ndarray
        Assembly Discontinuity Factors
    CDF : ndarray
        Corner Discontinuity Factors

    Raises
    ------
    RuntimeError
        If moc does not have a CMFD instance, or if it is not possible to
        create a condensation scheme to go from the CMFD group structure to
        the few-group structure.
    """
    if moc.cmfd is None:
        raise RuntimeError("moc must have a CMFD instance")

    if not moc.solved:
        raise RuntimeError("moc must be in a solved state")

    cmfd_condensation_scheme = moc.cmfd.condensation_scheme
    NG = len(moc_to_nodal_cond_scheme)

    # Create the condensation scheme to go from CMFD to few group structure
    cmfd_to_nodal_cond_scheme = []
    for G in range(len(moc_to_nodal_cond_scheme)):
        cmfd_grps = []

        # Get lower group
        for g in range(len(cmfd_condensation_scheme)):
            if cmfd_condensation_scheme[g][0] == moc_to_nodal_cond_scheme[G][0]:
                cmfd_grps.append(g)
        if len(cmfd_grps) == 0:
            raise RuntimeError(
                "Could not create condensation scheme from CMFD structure to few-group structure."
            )

        # Get upper group
        for g in range(len(cmfd_condensation_scheme)):
            if cmfd_condensation_scheme[g][-1] == moc_to_nodal_cond_scheme[G][-1]:
                cmfd_grps.append(g)
        if len(cmfd_grps) == 1:
            raise RuntimeError(
                "Could not create condensation scheme from CMFD structure to few-group structure."
            )

        cmfd_to_nodal_cond_scheme.append(cmfd_grps)

    # Create empty arrays for ADFs and CDFs
    adf = np.ones((NG, 6))
    cdf = np.ones((NG, 4))

    # Compute the homogeneous flux for each few-group
    hom_flux = _get_hom_flux_from_cmfd(
        moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
    )

    # Now we can go get the surface fluxes
    if symmetry in [Symmetry.Quarter, Symmetry.Half, Symmetry.Full]:
        het_flux_xp = _get_het_flux_xp_cmfd(
            moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
        )
        het_flux_yp = _get_het_flux_yp_cmfd(
            moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
        )
        het_flux_I = _get_het_flux_I_cmfd(
            moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
        )

        adf[:, ADF.XP] = het_flux_xp / hom_flux.flux_x.pos_surf_flux()
        adf[:, ADF.XN] = adf[:, ADF.XP]
        adf[:, ADF.YP] = het_flux_yp / hom_flux.flux_y.pos_surf_flux()
        adf[:, ADF.YN] = adf[:, ADF.YP]

        cdf[:, CDF.I] = het_flux_I / hom_flux(0.5 * hom_flux.dx, 0.5 * hom_flux.dy)
        cdf[:, CDF.II] = cdf[:, CDF.I]
        cdf[:, CDF.III] = cdf[:, CDF.I]
        cdf[:, CDF.IV] = cdf[:, CDF.I]

    if symmetry in [Symmetry.Half, Symmetry.Full] or independent_quadrant:
        het_flux_xn = _get_het_flux_xn_cmfd(
            moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
        )
        het_flux_II = _get_het_flux_II_cmfd(
            moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
        )

        adf[:, ADF.XN] = het_flux_xn / hom_flux.flux_x.neg_surf_flux()

        cdf[:, CDF.II] = het_flux_II / hom_flux(-0.5 * hom_flux.dx, 0.5 * hom_flux.dy)
        cdf[:, CDF.III] = cdf[:, CDF.II]

    if symmetry == Symmetry.Full or independent_quadrant:
        het_flux_yn = _get_het_flux_yn_cmfd(
            moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
        )
        het_flux_III = _get_het_flux_III_cmfd(
            moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
        )
        het_flux_IV = _get_het_flux_IV_cmfd(
            moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme
        )

        adf[:, ADF.YN] = het_flux_yn / hom_flux.flux_y.neg_surf_flux()

        cdf[:, CDF.III] = het_flux_III / hom_flux(
            -0.5 * hom_flux.dx, -0.5 * hom_flux.dy
        )
        cdf[:, CDF.IV] = het_flux_IV / hom_flux(0.5 * hom_flux.dx, -0.5 * hom_flux.dy)

    return adf, cdf


def compute_adf_cdf_from_moc(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    symmetry: Symmetry,
    independent_quadrant: bool,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the assembly and corner discontinuity factors in the few-group
    structure using the MOC results.

    Parameters
    ----------
    moc : MOCDriver
        Solved MOC problem for which the equivalence factors will be computed.
    moc_to_nodal_cond_scheme : list of list of ints
        Scheme used to condense from the MOC problem to the resulting nodal
        diffusion parameters.
    symmetry : Symmetry
        Indicates the symmetry of the problem.
    independent_quadrant : bool
        If quarter symmetry is indicated (symmetry = Symmetry.Quarter), then
        unique ADFs and CDFs will still be generated for each side, as if full
        symmetry were being used.

    Returns
    -------
    ADF : ndarray
        Assembly Discontinuity Factors
    CDF : ndarray
        Corner Discontinuity Factors
    """
    NG = len(moc_to_nodal_cond_scheme)

    # First, compute the homogeneous flux
    homog_flux = np.zeros(NG)
    total_volume = 0.0
    for i in range(moc.nfsr):
        Vi = moc.volume(i)
        total_volume += Vi

        for G in range(NG):
            gmin, gmax = moc_to_nodal_cond_scheme[G][:]
            for g in range(gmin, gmax + 1):
                homog_flux[G] += Vi * moc.flux(i, g)
    homog_flux /= total_volume

    # Create empty arrays for ADFs and CDFs
    adf = np.ones((NG, 6))
    cdf = np.zeros((NG, 4))

    if symmetry == Symmetry.Full or (
        symmetry == Symmetry.Quarter and independent_quadrant
    ):
        # Get flux along surfaces
        xn_segments = moc.trace_fsr_segments(
            Vector(moc.x_min + 0.001, moc.y_max), Direction(0.0, -1.0)
        )
        xp_segments = moc.trace_fsr_segments(
            Vector(moc.x_max - 0.001, moc.y_max), Direction(0.0, -1.0)
        )
        yn_segments = moc.trace_fsr_segments(
            Vector(moc.x_min, moc.y_min + 0.001), Direction(1.0, 0.0)
        )
        yp_segments = moc.trace_fsr_segments(
            Vector(moc.x_min, moc.y_max - 0.001), Direction(1.0, 0.0)
        )
        xn_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, xn_segments)
        xp_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, xp_segments)
        yn_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, yn_segments)
        yp_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, yp_segments)

        # Get flux at corners
        I_flx = _compute_few_group_flux(
            moc,
            moc_to_nodal_cond_scheme,
            Vector(moc.x_max - 0.001, moc.y_max - 0.001),
            Direction(-1.0, -1.0),
        )
        II_flx = _compute_few_group_flux(
            moc,
            moc_to_nodal_cond_scheme,
            Vector(moc.x_min + 0.001, moc.y_max - 0.001),
            Direction(1.0, -1.0),
        )
        III_flx = _compute_few_group_flux(
            moc,
            moc_to_nodal_cond_scheme,
            Vector(moc.x_min + 0.001, moc.y_min + 0.001),
            Direction(1.0, 1.0),
        )
        IV_flx = _compute_few_group_flux(
            moc,
            moc_to_nodal_cond_scheme,
            Vector(moc.x_max - 0.001, moc.y_min + 0.001),
            Direction(-1.0, 1.0),
        )

        for G in range(NG):
            adf[G, ADF.XN] = xn_flx[G] / homog_flux[G]
            adf[G, ADF.XP] = xp_flx[G] / homog_flux[G]
            adf[G, ADF.YN] = yn_flx[G] / homog_flux[G]
            adf[G, ADF.YP] = yp_flx[G] / homog_flux[G]
            # The ADFs on the +/- Z sides will be left at unity

            cdf[G, CDF.I] = I_flx[G] / homog_flux[G]
            cdf[G, CDF.II] = II_flx[G] / homog_flux[G]
            cdf[G, CDF.III] = III_flx[G] / homog_flux[G]
            cdf[G, CDF.IV] = IV_flx[G] / homog_flux[G]

    elif symmetry == Symmetry.Half:
        # Get flux along surfaces
        xn_segments = moc.trace_fsr_segments(
            Vector(moc.x_min + 0.001, moc.y_max), Direction(0.0, -1.0)
        )
        xp_segments = moc.trace_fsr_segments(
            Vector(moc.x_max - 0.001, moc.y_max), Direction(0.0, -1.0)
        )
        yp_segments = moc.trace_fsr_segments(
            Vector(moc.x_min, moc.y_max - 0.001), Direction(1.0, 0.0)
        )
        xn_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, xn_segments)
        xp_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, xp_segments)
        yp_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, yp_segments)

        # Get flux at corners
        I_flx = _compute_few_group_flux(
            moc,
            moc_to_nodal_cond_scheme,
            Vector(moc.x_max - 0.001, moc.y_max - 0.001),
            Direction(-1.0, -1.0),
        )
        II_flx = _compute_few_group_flux(
            moc,
            moc_to_nodal_cond_scheme,
            Vector(moc.x_min + 0.001, moc.y_max - 0.001),
            Direction(1.0, -1.0),
        )

        for G in range(NG):
            adf[G, ADF.XN] = xn_flx[G] / homog_flux[G]
            adf[G, ADF.XP] = xp_flx[G] / homog_flux[G]
            adf[G, ADF.YP] = yp_flx[G] / homog_flux[G]
            adf[G, ADF.YN] = adf[G, ADF.YP]
            # The ADFs on the +/- Z sides will be left at unity

            cdf[G, CDF.I] = I_flx[G] / homog_flux[G]
            cdf[G, CDF.II] = II_flx[G] / homog_flux[G]
            cdf[G, CDF.III] = cdf[G, CDF.II]
            cdf[G, CDF.IV] = cdf[G, CDF.I]

    else:  # Quarter symmetry
        # Get flux along surfaces
        xp_segments = moc.trace_fsr_segments(
            Vector(moc.x_max - 0.001, moc.y_max), Direction(0.0, -1.0)
        )
        yp_segments = moc.trace_fsr_segments(
            Vector(moc.x_min, moc.y_max - 0.001), Direction(1.0, 0.0)
        )
        xp_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, xp_segments)
        yp_flx = _compute_average_line_flux(moc, moc_to_nodal_cond_scheme, yp_segments)

        # Get flux at corners
        I_flx = _compute_few_group_flux(
            moc,
            moc_to_nodal_cond_scheme,
            Vector(moc.x_max - 0.001, moc.y_max - 0.001),
            Direction(-1.0, -1.0),
        )

        for G in range(NG):
            adf[G, ADF.XP] = xp_flx[G] / homog_flux[G]
            adf[G, ADF.XN] = adf[G, ADF.XP]
            adf[G, ADF.YP] = yp_flx[G] / homog_flux[G]
            adf[G, ADF.YN] = adf[G, ADF.YP]
            # The ADFs on the +/- Z sides will be left at unity

            cdf[G, CDF.I] = I_flx[G] / homog_flux[G]
            cdf[G, CDF.II] = cdf[G, CDF.I]
            cdf[G, CDF.III] = cdf[G, CDF.I]
            cdf[G, CDF.IV] = cdf[G, CDF.I]

    return adf, cdf


def _compute_few_group_flux(
    moc: MOCDriver, condensation_scheme: List[List[int]], r: Vector, u: Direction
) -> List[float]:
    """
    Computes the flux in the few-group scheme at the given position
    and direction.

    Parameters
    ----------
    moc : MOCDriver
        Solved MOC problem for which the flux will be computed.
    condensation_scheme : list of list of ints
        Scheme used to condense the flux estimates from the MOC problem.
    r : Vector
        Position where the flux is evaluated.
    u : Direction
        Direction vector to disambiguate the position.

    Returns
    -------
    list of float
        Values of the few-group flux at the given position.
    """
    flux = [0.0 for G in range(len(condensation_scheme))]

    for G in range(len(condensation_scheme)):
        gmin, gmax = condensation_scheme[G][:]

        for g in range(gmin, gmax + 1):
            flux[G] += moc.flux(r, u, g)

    return flux


def _compute_average_line_flux(
    moc: MOCDriver,
    condensation_scheme: List[List[int]],
    segments: List[Tuple[int, float]],
) -> np.ndarray:
    """
    Computes the average flux along a set of line segments in the few-group
    structure.

    Parameters
    ----------
    moc : MOCDriver
        Solved MOC problem for which the flux will be computed.
    condensation_scheme : list of list of ints
        Scheme used to condense the flux estimates from the MOC problem.
    segments : list of tuples of int and float
        List of flat source region index and segment length tuples.

    Returns
    -------
    ndarray
        Values of the few-group flux alone the line.
    """
    total_length = 0.0
    for s in segments:
        total_length += s[1]
    invs_tot_length = 1.0 / total_length

    flux = np.zeros(len(condensation_scheme))

    for G in range(len(condensation_scheme)):
        gmin, gmax = condensation_scheme[G][:]
        for g in range(gmin, gmax + 1):
            for s in segments:
                flux[G] += s[1] * moc.flux(s[0], g)

    flux *= invs_tot_length

    return flux


def _get_hom_flux_from_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
) -> NodalFlux2D:
    NG = len(cmfd_to_nodal_cond_scheme)
    keff = moc.keff

    dx = moc.x_max - moc.x_min
    dy = moc.y_max - moc.y_min

    # Compute average flux for assembly
    flux_spec = moc.homogenize_flux_spectrum()
    avg_flux = np.zeros(NG)
    for G in range(NG):
        for g in range(
            moc_to_nodal_cond_scheme[G][0], moc_to_nodal_cond_scheme[G][1] + 1
        ):
            avg_flux[G] += flux_spec[g]

    # Homogenize the xs and condense for assembly
    fine_xs = moc.homogenize()
    fine_diff_xs = fine_xs.diffusion_xs()
    few_diff_xs = fine_diff_xs.condense(moc_to_nodal_cond_scheme, flux_spec)

    # Condense currents for the assembly
    j_x_neg = _get_avg_x_neg_current_cmfd(moc, cmfd_to_nodal_cond_scheme)
    j_x_pos = _get_avg_x_pos_current_cmfd(moc, cmfd_to_nodal_cond_scheme)
    j_y_neg = _get_avg_y_neg_current_cmfd(moc, cmfd_to_nodal_cond_scheme)
    j_y_pos = _get_avg_y_pos_current_cmfd(moc, cmfd_to_nodal_cond_scheme)

    return NodalFlux2D(
        dx, dy, keff, few_diff_xs, avg_flux, j_x_neg, j_x_pos, j_y_neg, j_y_pos
    )


def _get_het_flux_I_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
) -> np.ndarray:
    cmfd = moc.cmfd
    i = cmfd.nx - 1
    j = cmfd.ny - 1

    tile_nodal_flux = _get_2d_nodal_flux_from_cmfd(
        moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme, i, j
    )
    return tile_nodal_flux(0.5 * tile_nodal_flux.dx, 0.5 * tile_nodal_flux.dy)


def _get_het_flux_II_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
) -> np.ndarray:
    cmfd = moc.cmfd
    i = 0
    j = cmfd.ny - 1

    tile_nodal_flux = _get_2d_nodal_flux_from_cmfd(
        moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme, i, j
    )
    return tile_nodal_flux(-0.5 * tile_nodal_flux.dx, 0.5 * tile_nodal_flux.dy)


def _get_het_flux_III_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
) -> np.ndarray:
    i = 0
    j = 0

    tile_nodal_flux = _get_2d_nodal_flux_from_cmfd(
        moc, moc_to_nodal_cond_scheme, cmfd_to_nodal_cond_scheme, i, j
    )
    return tile_nodal_flux(-0.5 * tile_nodal_flux.dx, -0.5 * tile_nodal_flux.dy)


def _get_het_flux_IV_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cond_scheme: List[List[int]],
) -> np.ndarray:
    cmfd = moc.cmfd
    i = cmfd.nx - 1
    j = 0

    tile_nodal_flux = _get_2d_nodal_flux_from_cmfd(
        moc, moc_to_nodal_cond_scheme, cond_scheme, i, j
    )
    return tile_nodal_flux(0.5 * tile_nodal_flux.dx, -0.5 * tile_nodal_flux.dy)


def _get_2d_nodal_flux_from_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
    i: int,
    j: int,
) -> NodalFlux2D:
    NG = len(moc_to_nodal_cond_scheme)
    keff = moc.keff
    cmfd = moc.cmfd

    tally_j_xp = i != cmfd.nx - 1 or moc.x_max_bc != BoundaryCondition.Reflective
    tally_j_xn = i != 0 or moc.x_min_bc != BoundaryCondition.Reflective
    tally_j_yp = j != cmfd.ny - 1 or moc.y_max_bc != BoundaryCondition.Reflective
    tally_j_yn = j != 0 or moc.y_min_bc != BoundaryCondition.Reflective

    dx = cmfd.dx[i]
    dy = cmfd.dy[j]

    # Get surface indices
    s_x_neg = cmfd.get_x_neg_surf(i, j)
    s_x_pos = cmfd.get_x_pos_surf(i, j)
    s_y_neg = cmfd.get_y_neg_surf(i, j)
    s_y_pos = cmfd.get_y_pos_surf(i, j)

    # Homogenize average flux for the tile
    cmfd_tile_fsrs = cmfd.tile_fsr_list(i, j)
    tile_flux_spec = moc.homogenize_flux_spectrum(cmfd_tile_fsrs)
    avg_flux = np.zeros(NG)
    for G in range(NG):
        for g in range(
            moc_to_nodal_cond_scheme[G][0], moc_to_nodal_cond_scheme[G][1] + 1
        ):
            avg_flux[G] += tile_flux_spec[g]

    # Homogenize the xs and condense
    fine_tile_xs = moc.homogenize(cmfd_tile_fsrs)
    fine_tile_diff_xs = fine_tile_xs.diffusion_xs()
    few_diff_xs = fine_tile_diff_xs.condense(moc_to_nodal_cond_scheme, tile_flux_spec)

    # Condense currents for the tile
    j_x_neg = np.zeros(NG)
    j_x_pos = np.zeros(NG)
    j_y_neg = np.zeros(NG)
    j_y_pos = np.zeros(NG)
    for G in range(NG):
        for g in range(
            cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
        ):
            if tally_j_xn:
                j_x_neg[G] += cmfd.current(g, s_x_neg)
            if tally_j_xp:
                j_x_pos[G] += cmfd.current(g, s_x_pos)
            if tally_j_yn:
                j_y_neg[G] += cmfd.current(g, s_y_neg)
            if tally_j_yp:
                j_y_pos[G] += cmfd.current(g, s_y_pos)

    return NodalFlux2D(
        dx, dy, keff, few_diff_xs, avg_flux, j_x_neg, j_x_pos, j_y_neg, j_y_pos
    )


def _get_avg_x_pos_current_cmfd(
    moc: MOCDriver, cmfd_to_nodal_cond_scheme: List[List[int]]
) -> np.ndarray:
    NG = len(cmfd_to_nodal_cond_scheme)
    cmfd = moc.cmfd

    current = np.zeros(NG)
    if moc.x_max_bc == BoundaryCondition.Reflective:
        return current

    dlts = cmfd.dy

    i = cmfd.nx - 1
    for j in range(cmfd.ny):
        # Get surface indices
        s = cmfd.get_x_pos_surf(i, j)

        # Condense currents for the tile
        for G in range(NG):
            for g in range(
                cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
            ):
                # Weight current contribution by length of segment
                current[G] += cmfd.current(g, s) * dlts[j]

    # Normalize by length of xp surface
    current /= moc.y_max - moc.y_min
    return current


def _get_avg_x_neg_current_cmfd(
    moc: MOCDriver, cmfd_to_nodal_cond_scheme: List[List[int]]
) -> np.ndarray:
    NG = len(cmfd_to_nodal_cond_scheme)
    cmfd = moc.cmfd

    current = np.zeros(NG)
    if moc.x_min_bc == BoundaryCondition.Reflective:
        return current

    dlts = cmfd.dy

    i = 0
    for j in range(cmfd.ny):
        # Get surface indices
        s = cmfd.get_x_neg_surf(i, j)

        # Condense currents for the tile
        for G in range(NG):
            for g in range(
                cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
            ):
                # Weight current contribution by length of segment
                current[G] += cmfd.current(g, s) * dlts[j]

    # Normalize by length of xp surface
    current /= moc.y_max - moc.y_min
    return current


def _get_avg_y_pos_current_cmfd(
    moc: MOCDriver, cmfd_to_nodal_cond_scheme: List[List[int]]
) -> np.ndarray:
    NG = len(cmfd_to_nodal_cond_scheme)
    cmfd = moc.cmfd

    current = np.zeros(NG)
    if moc.y_max_bc == BoundaryCondition.Reflective:
        return current

    dlts = cmfd.dx

    j = cmfd.ny - 1
    for i in range(cmfd.nx):
        # Get surface indices
        s = cmfd.get_y_pos_surf(i, j)

        # Condense currents for the tile
        for G in range(NG):
            for g in range(
                cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
            ):
                # Weight current contribution by length of segment
                current[G] += cmfd.current(g, s) * dlts[i]

    # Normalize by length of xp surface
    current /= moc.x_max - moc.x_min
    return current


def _get_avg_y_neg_current_cmfd(
    moc: MOCDriver, cmfd_to_nodal_cond_scheme: List[List[int]]
) -> np.ndarray:
    NG = len(cmfd_to_nodal_cond_scheme)
    cmfd = moc.cmfd

    current = np.zeros(NG)
    if moc.y_min_bc == BoundaryCondition.Reflective:
        return current

    dlts = cmfd.dx

    j = 0
    for i in range(cmfd.nx):
        # Get surface indices
        s = cmfd.get_y_neg_surf(i, j)

        # Condense currents for the tile
        for G in range(NG):
            for g in range(
                cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
            ):
                # Weight current contribution by length of segment
                current[G] += cmfd.current(g, s) * dlts[i]

    # Normalize by length of xp surface
    current /= moc.x_max - moc.x_min
    return current


def _get_het_flux_xp_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
) -> np.ndarray:
    NG = len(cmfd_to_nodal_cond_scheme)
    keff = moc.keff
    cmfd = moc.cmfd

    tally_j_pos = moc.x_max_bc != BoundaryCondition.Reflective

    dlts = cmfd.dy

    node_fluxes = []

    i = cmfd.nx - 1
    x_min = 0.0
    x_max = cmfd.dx[-1]
    for j in range(cmfd.ny):
        # For tile (i, j), get the FSR list
        cmfd_tile_fsrs = cmfd.tile_fsr_list(i, j)

        # Get surface indices
        s_neg = cmfd.get_x_neg_surf(i, j)
        s_pos = cmfd.get_x_pos_surf(i, j)

        # Homogenize average flux for the tile
        tile_flux_spec = moc.homogenize_flux_spectrum(cmfd_tile_fsrs)
        avg_flux = np.zeros(NG)
        for G in range(NG):
            for g in range(
                moc_to_nodal_cond_scheme[G][0], moc_to_nodal_cond_scheme[G][1] + 1
            ):
                avg_flux[G] += tile_flux_spec[g]

        # Homogenize the xs and condense
        fine_tile_xs = moc.homogenize(cmfd_tile_fsrs)
        fine_tile_diff_xs = fine_tile_xs.diffusion_xs()
        few_diff_xs = fine_tile_diff_xs.condense(
            moc_to_nodal_cond_scheme, tile_flux_spec
        )

        # Condense currents for the tile
        j_neg = np.zeros(NG)
        j_pos = np.zeros(NG)
        for G in range(NG):
            for g in range(
                cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
            ):
                j_neg[G] += cmfd.current(g, s_neg)
                if tally_j_pos:
                    j_pos[G] += cmfd.current(g, s_pos)

        # Create 1D nodal flux object
        node_fluxes.append(
            NodalFlux1D(x_min, x_max, keff, few_diff_xs, avg_flux, j_neg, j_pos)
        )

    # Now we need to get the average heterogeneous flux on the surface
    het_flux = np.zeros(NG)
    for G in range(NG):
        for j in range(cmfd.ny):
            het_flux[G] += node_fluxes[j].pos_surf_flux(G) * dlts[j]
    het_flux /= moc.y_max - moc.y_min

    return het_flux


def _get_het_flux_xn_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
) -> np.ndarray:
    NG = len(cmfd_to_nodal_cond_scheme)
    keff = moc.keff
    cmfd = moc.cmfd

    tally_j_neg = moc.x_min_bc != BoundaryCondition.Reflective

    dlts = cmfd.dy

    node_fluxes = []

    i = 0
    x_min = 0.0
    x_max = cmfd.dx[0]
    for j in range(cmfd.ny):
        # For tile (i, j), get the FSR list
        cmfd_tile_fsrs = cmfd.tile_fsr_list(i, j)

        # Get surface indices
        s_neg = cmfd.get_x_neg_surf(i, j)
        s_pos = cmfd.get_x_pos_surf(i, j)

        # Homogenize average flux for the tile
        tile_flux_spec = moc.homogenize_flux_spectrum(cmfd_tile_fsrs)
        avg_flux = np.zeros(NG)
        for G in range(NG):
            for g in range(
                moc_to_nodal_cond_scheme[G][0], moc_to_nodal_cond_scheme[G][1] + 1
            ):
                avg_flux[G] += tile_flux_spec[g]

        # Homogenize the xs and condense
        fine_tile_xs = moc.homogenize(cmfd_tile_fsrs)
        fine_tile_diff_xs = fine_tile_xs.diffusion_xs()
        few_diff_xs = fine_tile_diff_xs.condense(
            moc_to_nodal_cond_scheme, tile_flux_spec
        )

        # Condense currents for the tile
        j_neg = np.zeros(NG)
        j_pos = np.zeros(NG)
        for G in range(NG):
            for g in range(
                cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
            ):
                if tally_j_neg:
                    j_neg[G] += cmfd.current(g, s_neg)
                j_pos[G] += cmfd.current(g, s_pos)

        # Create 1D nodal flux object
        node_fluxes.append(
            NodalFlux1D(x_min, x_max, keff, few_diff_xs, avg_flux, j_neg, j_pos)
        )

    # Now we need to get the average heterogeneous flux on the surface
    het_flux = np.zeros(NG)
    for G in range(NG):
        for j in range(cmfd.ny):
            het_flux[G] += node_fluxes[j].neg_surf_flux(G) * dlts[j]
    het_flux /= moc.y_max - moc.y_min

    return het_flux


def _get_het_flux_yp_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
) -> np.ndarray:
    NG = len(cmfd_to_nodal_cond_scheme)
    keff = moc.keff
    cmfd = moc.cmfd

    tally_j_pos = moc.y_max_bc != BoundaryCondition.Reflective

    dlts = cmfd.dx

    node_fluxes = []

    j = cmfd.ny - 1
    y_min = 0.0
    y_max = cmfd.dy[-1]
    for i in range(cmfd.nx):
        # For tile (i, j), get the FSR list
        cmfd_tile_fsrs = cmfd.tile_fsr_list(i, j)

        # Get surface indices
        s_neg = cmfd.get_y_neg_surf(i, j)
        s_pos = cmfd.get_y_pos_surf(i, j)

        # Homogenize average flux for the tile
        tile_flux_spec = moc.homogenize_flux_spectrum(cmfd_tile_fsrs)
        avg_flux = np.zeros(NG)
        for G in range(NG):
            for g in range(
                moc_to_nodal_cond_scheme[G][0], moc_to_nodal_cond_scheme[G][1] + 1
            ):
                avg_flux[G] += tile_flux_spec[g]

        # Homogenize the xs and condense
        fine_tile_xs = moc.homogenize(cmfd_tile_fsrs)
        fine_tile_diff_xs = fine_tile_xs.diffusion_xs()
        few_diff_xs = fine_tile_diff_xs.condense(
            moc_to_nodal_cond_scheme, tile_flux_spec
        )

        # Condense currents for the tile
        j_neg = np.zeros(NG)
        j_pos = np.zeros(NG)
        for G in range(NG):
            for g in range(
                cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
            ):
                j_neg[G] += cmfd.current(g, s_neg)
                if tally_j_pos:
                    j_pos[G] += cmfd.current(g, s_pos)

        # Create 1D nodal flux object
        node_fluxes.append(
            NodalFlux1D(y_min, y_max, keff, few_diff_xs, avg_flux, j_neg, j_pos)
        )

    # Now we need to get the average heterogeneous flux on the surface
    het_flux = np.zeros(NG)
    for G in range(NG):
        for i in range(cmfd.nx):
            het_flux[G] += node_fluxes[i].pos_surf_flux(G) * dlts[i]
    het_flux /= moc.x_max - moc.x_min

    return het_flux


def _get_het_flux_yn_cmfd(
    moc: MOCDriver,
    moc_to_nodal_cond_scheme: List[List[int]],
    cmfd_to_nodal_cond_scheme: List[List[int]],
) -> np.ndarray:
    NG = len(cmfd_to_nodal_cond_scheme)
    keff = moc.keff
    cmfd = moc.cmfd

    tally_j_neg = moc.y_min_bc != BoundaryCondition.Reflective

    dlts = cmfd.dx
    node_fluxes = []

    j = 0
    y_min = 0.0
    y_max = cmfd.dy[0]
    for i in range(cmfd.nx):
        # For tile (i, j), get the FSR list
        cmfd_tile_fsrs = cmfd.tile_fsr_list(i, j)

        # Get surface indices
        s_neg = cmfd.get_y_neg_surf(i, j)
        s_pos = cmfd.get_y_pos_surf(i, j)

        # Homogenize average flux for the tile
        tile_flux_spec = moc.homogenize_flux_spectrum(cmfd_tile_fsrs)
        avg_flux = np.zeros(NG)
        for G in range(NG):
            for g in range(
                moc_to_nodal_cond_scheme[G][0], moc_to_nodal_cond_scheme[G][1] + 1
            ):
                avg_flux[G] += tile_flux_spec[g]

        # Homogenize the xs and condense
        fine_tile_xs = moc.homogenize(cmfd_tile_fsrs)
        fine_tile_diff_xs = fine_tile_xs.diffusion_xs()
        few_diff_xs = fine_tile_diff_xs.condense(
            moc_to_nodal_cond_scheme, tile_flux_spec
        )

        # Condense currents for the tile
        j_neg = np.zeros(NG)
        j_pos = np.zeros(NG)
        for G in range(NG):
            for g in range(
                cmfd_to_nodal_cond_scheme[G][0], cmfd_to_nodal_cond_scheme[G][1] + 1
            ):
                if tally_j_neg:
                    j_neg[G] += cmfd.current(g, s_neg)
                j_pos[G] += cmfd.current(g, s_pos)

        # Create 1D nodal flux object
        node_fluxes.append(
            NodalFlux1D(y_min, y_max, keff, few_diff_xs, avg_flux, j_neg, j_pos)
        )

    # Now we need to get the average heterogeneous flux on the surface
    het_flux = np.zeros(NG)
    for G in range(NG):
        for i in range(cmfd.nx):
            het_flux[G] += node_fluxes[i].neg_surf_flux(G) * dlts[i]
    het_flux /= moc.x_max - moc.x_min

    return het_flux
