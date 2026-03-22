from pathlib import Path
import numpy as np
import libraries
from libraries import Nuclide, ThermalScatteringMaterial
import tools.frendy as fdy
import tools.depletion_chain as dc
import h5py
import logging
from scarabee import nuclide_name_to_za, nuclide_name_to_element_symbol

# Set up logging to file and console
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.FileHandler("generate_nuclear_data.log"),
                        logging.StreamHandler()
                    ])


####################################################################################
##                         MODIFY LINES IN THIS BLOCK

library = libraries.endf80_library # Change to desired library
library.download_library_if_not_exists()
# or, define custom library:
# library = libraries.ENDFLibrary(
#     name="Custom Library",
#     label="custom_label", # one of "endf71", "endf80", "endf81", "jeff33", "jeff40"
#     base_path=Path("/path/to/library"),
#     neutrons_path_suffix="neutrons",
#     tsl_path_suffix="thermal_scatt",
# )

temps = [293.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0]

hdf5_file_name = f"{library.label}_shem281.h5"

# Set the default group strucutre
fdy.set_default_group_structure("SHEM-281")
fdy.set_default_max_legendre_moments(3)

####################################################################################


# First we initialize the file with the depletion chain data.
# Must be done first as the save method will delete file first if it already exists !
dep_chain = dc.build_depletion_chain("chain_casl_pwr.xml")
dep_chain.save(hdf5_file_name)

# Initialize HDF5 with the library information
print(fdy.get_default_group_structure())
h5 = h5py.File(hdf5_file_name, "r+")
h5.attrs["group-structure"] = fdy.get_default_group_structure().name
h5.attrs["group-bounds"] = fdy.get_default_group_structure().bounds
h5.attrs["ngroups"] = fdy.get_default_group_structure().ngroups
h5.attrs["condensation-scheme"] = fdy.get_default_group_structure().condensation_scheme
h5.attrs["cmfd-condensation-scheme"] = (
    fdy.get_default_group_structure().cmfd_condensation_scheme
)
h5.attrs["first-resonance-group"] = fdy.get_default_group_structure().first_res_grp
h5.attrs["last-resonance-group"] = fdy.get_default_group_structure().last_res_grp
h5.attrs["library"] = library.name

# Define a helper function to process nuclides
def process_nuclide(nuclide: Nuclide, dilutions: list[float] = [1.0e10]):
    endf_file_path = library.neutrons_path / library.get_filename_of_nuclide(nuclide)
    if not endf_file_path.exists():
        logging.warning(f"ENDF file for {nuclide.label()} not found in library (looked for it at '{endf_file_path}')")
        return

    N = fdy.FrendyMG()
    N.name = nuclide.label()
    N.endf_file = str(endf_file_path)
    N.label = f"{nuclide.label()} from {library.name}"
    N.temps = temps
    if dilutions:
        N.dilutions = dilutions
    N.process(h5)

    logging.info(f"Processed nuclide: {nuclide.label()}")

    return N

def convert_nuclide_name_to_nuclide(nuclide_name: str, is_fission_product: bool = False) -> Nuclide:
    # nuclide_name_to_za returns (Z * 1000 + A) * 10 + m;
    za = nuclide_name_to_za(nuclide_name)

    Z = za // 10000
    A = (za % 10000) // 10
    m = za % 10

    if m == 0:
        m = None

    element_symbol = nuclide_name_to_element_symbol(nuclide_name)

    return Nuclide(
        symbol=element_symbol,
        Z=Z,
        A=A,
        is_fission_product=is_fission_product,
        metastable_state=m
    )

# Define a helper function to convert nuclide names to Nuclide objects
# with a shorter syntax for use in the list of nuclides to process
parse_nuclide = lambda name, is_fission_product=False: convert_nuclide_name_to_nuclide(name, is_fission_product=is_fission_product)

# Define a helper function to process TSL materials
def process_tsl_material(tsl_material: ThermalScatteringMaterial, nuclide: Nuclide, chi: np.ndarray | None, dilutions: list[float] = [1.0e10]):
    tsl_temps = library.get_tsl_temperatures(tsl_material)
    
    endf_file_path = library.neutrons_path / library.get_filename_of_nuclide(nuclide)
    if not endf_file_path.exists():
        logging.warning(f"ENDF file for {nuclide.label()} not found in library (looked for it at '{endf_file_path}')")
        return
    tsl_endf_file_path = library.tsl_path / library.get_filename_of_tsl_material(tsl_material)
    if not tsl_endf_file_path.exists():
        logging.warning(f"ENDF TSL file for {tsl_material.label()} not found in library (looked for it at '{tsl_endf_file_path}')")
        return

    N = fdy.FrendyMG()
    N.name = tsl_material.label()
    N.endf_file = str(endf_file_path)
    N.tsl_file = str(tsl_endf_file_path)
    N.tsl_type = tsl_material.type.lower()
    N.label = f"{tsl_material.label()} from {library.name}"
    N.temps = tsl_temps
    if dilutions:
        N.dilutions = dilutions
    N.process(h5, chi)

    logging.info(f"Processed TSL material: {tsl_material.label()} for nuclide {nuclide.label()}")

    return N

# We process U235 first so that we can steal its fission spectrum for
# performing transport correciton calculations
N = process_nuclide(parse_nuclide("U235"))

if not N:
    raise RuntimeError("U235 data could not be processed")

if N.chi is not None:
    chi = N.chi[:]
else:
    raise RuntimeError("U235 fission spectrum not found")

# Process TSL based evaluations
process_tsl_material(
    ThermalScatteringMaterial("H", "H2O", "hh2o"),
    parse_nuclide("H1"),
    chi=chi,
    dilutions=[1.0e10],
)

process_tsl_material(
    ThermalScatteringMaterial("D", "D2O", "dd2o"),
    parse_nuclide("H2"),
    chi=chi,
    dilutions=[1.0e10],
)

# Free Gas Nuclides
nuclides = [
    parse_nuclide("H1"),
    parse_nuclide("H2"),
    parse_nuclide("He3"),
    parse_nuclide("He4"),
    parse_nuclide("Li6"),
    parse_nuclide("Li7"),
    parse_nuclide("Be9"),
    parse_nuclide("B10"),
    parse_nuclide("B11"),
    parse_nuclide("C12"),
    parse_nuclide("C13"),
    parse_nuclide("N14"),
    parse_nuclide("N15"),
    parse_nuclide("O16"),
    parse_nuclide("O17"),
    parse_nuclide("O18"),
    parse_nuclide("Na23"),
    parse_nuclide("Mg24"),
    parse_nuclide("Mg25"),
    parse_nuclide("Mg26"),
    parse_nuclide("Al27"),
    parse_nuclide("Si28"),
    parse_nuclide("Si29"),
    parse_nuclide("Si30"),
    parse_nuclide("Ar36"),
    parse_nuclide("Ar38"),
    parse_nuclide("Ar40"),
    parse_nuclide("Ti46"),
    parse_nuclide("Ti47"),
    parse_nuclide("Ti48"),
    parse_nuclide("Ti49"),
    parse_nuclide("Ti50"),
    parse_nuclide("Cr50"),
    parse_nuclide("Cr52"),
    parse_nuclide("Cr53"),
    parse_nuclide("Cr54"),
    parse_nuclide("Mn55"),
    parse_nuclide("Fe54"),
    parse_nuclide("Fe55"),
    parse_nuclide("Fe56"),
    parse_nuclide("Fe57"),
    parse_nuclide("Fe58"),
    parse_nuclide("Co59"),
    parse_nuclide("Ni58"),
    parse_nuclide("Ni60"),
    parse_nuclide("Ni61"),
    parse_nuclide("Ni62"),
    parse_nuclide("Ni64"),
    parse_nuclide("Cu63"),
    parse_nuclide("Cu65"),
    parse_nuclide("Zn64"),
    parse_nuclide("Zn66"),
    parse_nuclide("Zn67"),
    parse_nuclide("Zn68"),
    parse_nuclide("Zn70"),
    parse_nuclide("Br81"),
    parse_nuclide("Kr82"),
    parse_nuclide("Kr83"),
    parse_nuclide("Kr84"),
    parse_nuclide("Kr85"),
    parse_nuclide("Kr86"),
    parse_nuclide("Sr89"),
    parse_nuclide("Sr90"),
    parse_nuclide("Y89"),
    parse_nuclide("Y90"),
    parse_nuclide("Y91"),
    parse_nuclide("Zr90"),
    parse_nuclide("Zr91"),
    parse_nuclide("Zr92"),
    parse_nuclide("Zr93", is_fission_product=True),
    parse_nuclide("Zr94"),
    parse_nuclide("Zr95", is_fission_product=True),
    parse_nuclide("Zr96"),
    parse_nuclide("Nb95", is_fission_product=True),
    parse_nuclide("Mo92"),
    parse_nuclide("Mo94"),
    parse_nuclide("Mo95"),
    parse_nuclide("Mo96"),
    parse_nuclide("Mo97"),
    parse_nuclide("Mo98"),
    parse_nuclide("Mo99", is_fission_product=True),
    parse_nuclide("Mo100"),
    parse_nuclide("Tc99"),
    parse_nuclide("Ru99", is_fission_product=True),
    parse_nuclide("Ru100", is_fission_product=True),
    parse_nuclide("Ru101", is_fission_product=True),
    parse_nuclide("Ru102", is_fission_product=True),
    parse_nuclide("Ru103", is_fission_product=True),
    parse_nuclide("Ru104", is_fission_product=True),
    parse_nuclide("Ru105", is_fission_product=True),
    parse_nuclide("Ru106", is_fission_product=True),
    parse_nuclide("Rh103", is_fission_product=True),
    parse_nuclide("Rh105", is_fission_product=True),
    parse_nuclide("Pd104", is_fission_product=True),
    parse_nuclide("Pd105", is_fission_product=True),
    parse_nuclide("Pd106", is_fission_product=True),
    parse_nuclide("Pd107", is_fission_product=True),
    parse_nuclide("Pd108", is_fission_product=True),
    parse_nuclide("Ag107"),
    parse_nuclide("Ag109"),
    parse_nuclide("Ag110m1"),
    parse_nuclide("Ag111"),
    parse_nuclide("Cd106"),
    parse_nuclide("Cd108"),
    parse_nuclide("Cd110"),
    parse_nuclide("Cd111"),
    parse_nuclide("Cd112"),
    parse_nuclide("Cd113"),
    parse_nuclide("Cd114"),
    parse_nuclide("Cd116"),
    parse_nuclide("In113"),
    parse_nuclide("In115"),
    parse_nuclide("Sn112"),
    parse_nuclide("Sn114"),
    parse_nuclide("Sn115"),
    parse_nuclide("Sn116"),
    parse_nuclide("Sn117"),
    parse_nuclide("Sn118"),
    parse_nuclide("Sn119"),
    parse_nuclide("Sn120"),
    parse_nuclide("Sn122"),
    parse_nuclide("Sn124"),
    parse_nuclide("Sb121", is_fission_product=True),
    parse_nuclide("Sb123", is_fission_product=True),
    parse_nuclide("Sb125", is_fission_product=True),
    parse_nuclide("Te127m1", is_fission_product=True),
    parse_nuclide("Te129m1", is_fission_product=True),
    parse_nuclide("Te132", is_fission_product=True),
    parse_nuclide("I127"),
    parse_nuclide("I129"),
    parse_nuclide("I131", is_fission_product=True),
    parse_nuclide("I135", is_fission_product=True),
    parse_nuclide("Xe128", is_fission_product=True),
    parse_nuclide("Xe129"),
    parse_nuclide("Xe130", is_fission_product=True),
    parse_nuclide("Xe131"),
    parse_nuclide("Xe132", is_fission_product=True),
    parse_nuclide("Xe133", is_fission_product=True),
    parse_nuclide("Xe134", is_fission_product=True),
    parse_nuclide("Xe135", is_fission_product=True),
    parse_nuclide("Xe136", is_fission_product=True),
    parse_nuclide("Cs133", is_fission_product=True),
    parse_nuclide("Cs134", is_fission_product=True),
    parse_nuclide("Cs135", is_fission_product=True),
    parse_nuclide("Cs136", is_fission_product=True),
    parse_nuclide("Cs137", is_fission_product=True),
    parse_nuclide("Ba134", is_fission_product=True),
    parse_nuclide("Ba137", is_fission_product=True),
    parse_nuclide("Ba140", is_fission_product=True),
    parse_nuclide("La139", is_fission_product=True),
    parse_nuclide("La140", is_fission_product=True),
    parse_nuclide("Ce140", is_fission_product=True),
    parse_nuclide("Ce141", is_fission_product=True),
    parse_nuclide("Ce142", is_fission_product=True),
    parse_nuclide("Ce143", is_fission_product=True),
    parse_nuclide("Ce144", is_fission_product=True),
    parse_nuclide("Pr141", is_fission_product=True),
    parse_nuclide("Pr143", is_fission_product=True),
    parse_nuclide("Nd142", is_fission_product=True),
    parse_nuclide("Nd143", is_fission_product=True),
    parse_nuclide("Nd144", is_fission_product=True),
    parse_nuclide("Nd145", is_fission_product=True),
    parse_nuclide("Nd146", is_fission_product=True),
    parse_nuclide("Nd147", is_fission_product=True),
    parse_nuclide("Nd148", is_fission_product=True),
    parse_nuclide("Nd150", is_fission_product=True),
    parse_nuclide("Pm147", is_fission_product=True),
    parse_nuclide("Pm148", is_fission_product=True),
    parse_nuclide("Pm148m1", is_fission_product=True),
    parse_nuclide("Pm149", is_fission_product=True),
    parse_nuclide("Pm151", is_fission_product=True),
    parse_nuclide("Sm147"),
    parse_nuclide("Sm148", is_fission_product=True),
    parse_nuclide("Sm149"),
    parse_nuclide("Sm150"),
    parse_nuclide("Sm151"),
    parse_nuclide("Sm152"),
    parse_nuclide("Sm153"),
    parse_nuclide("Sm154"),
    parse_nuclide("Eu151"),
    parse_nuclide("Eu152"),
    parse_nuclide("Eu153"),
    parse_nuclide("Eu154"),
    parse_nuclide("Eu155"),
    parse_nuclide("Eu156", is_fission_product=True),
    parse_nuclide("Gd152"),
    parse_nuclide("Gd154"),
    parse_nuclide("Gd155"),
    parse_nuclide("Gd156"),
    parse_nuclide("Gd157"),
    parse_nuclide("Gd158"),
    parse_nuclide("Gd160"),
    parse_nuclide("Tb159", is_fission_product=True),
    parse_nuclide("Tb160", is_fission_product=True),
    parse_nuclide("Tb161", is_fission_product=True),
    parse_nuclide("Dy160"),
    parse_nuclide("Dy161"),
    parse_nuclide("Dy162"),
    parse_nuclide("Dy163"),
    parse_nuclide("Dy164"),
    parse_nuclide("Ho165", is_fission_product=True),
    parse_nuclide("Er162", is_fission_product=True),
    parse_nuclide("Er164", is_fission_product=True),
    parse_nuclide("Er166", is_fission_product=True),
    parse_nuclide("Er167", is_fission_product=True),
    parse_nuclide("Er168", is_fission_product=True),
    parse_nuclide("Er169", is_fission_product=True),
    parse_nuclide("Er170", is_fission_product=True),
    parse_nuclide("Tm169", is_fission_product=True),
    parse_nuclide("Tm170", is_fission_product=True),
    parse_nuclide("Tm171", is_fission_product=True),
    parse_nuclide("Hf174"),
    parse_nuclide("Hf176"),
    parse_nuclide("Hf177"),
    parse_nuclide("Hf178"),
    parse_nuclide("Hf179"),
    parse_nuclide("Hf180"),
    parse_nuclide("Hf181"),
    parse_nuclide("Ta181"),
    parse_nuclide("Ta182", is_fission_product=True),
    parse_nuclide("Th230"),
    parse_nuclide("Th231"),
    parse_nuclide("Th232"),
    parse_nuclide("Th234"),
    parse_nuclide("Pa231"),
    parse_nuclide("Pa232"),
    parse_nuclide("Pa233"),
    parse_nuclide("U232"),
    parse_nuclide("U233"),
    parse_nuclide("U234"),
    parse_nuclide("U236"),
    parse_nuclide("U237"),
    parse_nuclide("U238"),
    parse_nuclide("Np236"),
    parse_nuclide("Np237"),
    parse_nuclide("Np238"),
    parse_nuclide("Np239"),
    parse_nuclide("Pu236"),
    parse_nuclide("Pu237"),
    parse_nuclide("Pu238"),
    parse_nuclide("Pu239"),
    parse_nuclide("Pu240"),
    parse_nuclide("Pu241"),
    parse_nuclide("Pu242"),
    parse_nuclide("Pu243"),
    parse_nuclide("Pu244"),
    parse_nuclide("Am241"),
    parse_nuclide("Am242m1"),
    parse_nuclide("Am243"),
    parse_nuclide("Cm242"),
    parse_nuclide("Cm243"),
    parse_nuclide("Cm244"),
    parse_nuclide("Cm245"),
    parse_nuclide("Cm246"),
]

# Handle elemental evaluations
# See https://github.com/scarabee-dev/scarabee/pull/19#discussion_r2616262546
match library.label:
    case "endf71":
        # Remove all C nuclides and add elemental C evaluation
        nuclides = [nuclide for nuclide in nuclides if nuclide.symbol != "C"]
        nuclides.append(parse_nuclide("C0"))
    case "jeff33":
        # Remove C12 evaluation and add elemental C evaluation
        nuclides = [nuclide for nuclide in nuclides if not (nuclide.symbol == "C" and nuclide.A == 12)]
        nuclides.append(parse_nuclide("C0"))

# TODO: Add parallelization
[process_nuclide(nuclide) for nuclide in nuclides]

h5.close()
