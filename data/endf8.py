import tools.frendy as fdy
import tools.depletion_chain as dc
import h5py


####################################################################################
##                         MODIFY LINES IN THIS BLOCK
base = "/path/to/ENDF-B-VIII.0/neutrons/"
tslbase = "/path/to/ENDF-B-VIII.0/thermal_scatt/"
lib_name = "ENDF/B-VIII.0"

temps = [293.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0]

hdf5_file_name = "endf8_shem281.h5"

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
h5.attrs["library"] = lib_name

# Define a helper function to process nuclides
def process_nuclide(name: str, endf_filename: str, dilutions: list[float] | None = None):
    N = fdy.FrendyMG()
    N.name = name
    N.endf_file = base + endf_filename
    N.label = N.name + " from ENDF/B-8.0"
    N.temps = temps
    if dilutions:
        N.dilutions = dilutions
    N.process(h5)
    return N

# We process U235 first so that we can steal its fission spectrum for
# performing transport correciton calculations
N = process_nuclide("U235", "n-092_U_235.endf")

if N.chi:
    chi = N.chi[:]
else:
    raise RuntimeError("U235 fission spectrum not found")

# Process TSL based evaluations
N = fdy.FrendyMG()
N.name = "H1_H2O"
N.endf_file = base + "n-001_H_001.endf"
N.tsl_file = tslbase + "tsl-HinH2O.endf"
N.tsl_type = "hh2o"
N.label = "H1 in H2O from ENDF/B-8.0"
N.temps = [
    283.6,
    293.6,
    300.0,
    323.6,
    350.0,
    373.6,
    400.0,
    423.6,
    450.0,
    473.6,
    500.0,
    523.6,
    550.0,
    573.6,
    600.0,
    623.6,
    650.0,
    800.0,
]
N.dilutions = [1.0e10]
N.process(h5, chi)  # Processed with in-scatter transport correction !

N = fdy.FrendyMG()
N.name = "H2_D2O"
N.endf_file = base + "n-001_H_002.endf"
N.tsl_file = tslbase + "tsl-DinD2O.endf"
N.tsl_type = "dd2o"
N.label = "H2 in D2O from ENDF/B-8.0"
N.temps = [
    283.6,
    293.6,
    300.0,
    323.6,
    350.0,
    373.6,
    400.0,
    423.6,
    450.0,
    473.6,
    500.0,
    523.6,
    550.0,
    573.6,
    600.0,
    623.6,
    650.0,
]
N.dilutions = [1.0e10]
N.process(h5, chi)  # Processed with in-scatter transport correction !

# Free Gas Nuclides
process_nuclide("H1", "n-001_H_001.endf", [1.0e10])
process_nuclide("H1", "n-001_H_001.endf", [1.0e10])
process_nuclide("H2", "n-001_H_002.endf", [1.0e10])
process_nuclide("He3", "n-002_He_003.endf", [1.0e10])
process_nuclide("He4", "n-002_He_004.endf", [1.0e10])
process_nuclide("Li6", "n-003_Li_006.endf", [1.0e10])
process_nuclide("Li7", "n-003_Li_007.endf", [1.0e10])
process_nuclide("Be9", "n-004_Be_009.endf", [1.0e10])
process_nuclide("B10", "n-005_B_010.endf", [1.0e10])
process_nuclide("B11", "n-005_B_011.endf", [1.0e10])
process_nuclide("C12", "n-006_C_012.endf", [1.0e10])
process_nuclide("C13", "n-006_C_013.endf", [1.0e10])
process_nuclide("N14", "n-007_N_014.endf", [1.0e10])
process_nuclide("N15", "n-007_N_015.endf", [1.0e10])
process_nuclide("O16", "n-008_O_016.endf", [1.0e10])
process_nuclide("O17", "n-008_O_017.endf", [1.0e10])
process_nuclide("O18", "n-008_O_018.endf", [1.0e10])
process_nuclide("Na23", "n-011_Na_023.endf", [1.0e10])
process_nuclide("Mg24", "n-012_Mg_024.endf", [1.0e10])
process_nuclide("Mg25", "n-012_Mg_025.endf", [1.0e10])
process_nuclide("Mg26", "n-012_Mg_026.endf", [1.0e10])
process_nuclide("Al27", "n-013_Al_027.endf", [1.0e10])
process_nuclide("Si28", "n-014_Si_028.endf", [1.0e10])
process_nuclide("Si29", "n-014_Si_029.endf", [1.0e10])
process_nuclide("Si30", "n-014_Si_030.endf", [1.0e10])
process_nuclide("Ar36", "n-018_Ar_036.endf", [1.0e10])
process_nuclide("Ar38", "n-018_Ar_038.endf", [1.0e10])
process_nuclide("Ar40", "n-018_Ar_040.endf", [1.0e10])
process_nuclide("Cr50", "n-024_Cr_050.endf", [1.0e10])
process_nuclide("Cr52", "n-024_Cr_052.endf")
process_nuclide("Cr53", "n-024_Cr_053.endf")
process_nuclide("Cr54", "n-024_Cr_054.endf", [1.0e10])
process_nuclide("Mn55", "n-025_Mn_055.endf", [1.0e10])
process_nuclide("Fe54", "n-026_Fe_054.endf")
process_nuclide("Fe55", "n-026_Fe_055.endf", [1.0e10])
process_nuclide("Fe56", "n-026_Fe_056.endf")
process_nuclide("Fe57", "n-026_Fe_057.endf", [1.0e10])
process_nuclide("Fe58", "n-026_Fe_058.endf", [1.0e10])
process_nuclide("Co59", "n-027_Co_059.endf", [1.0e10])
process_nuclide("Ni58", "n-028_Ni_058.endf", [1.0e10])
process_nuclide("Ni60", "n-028_Ni_060.endf", [1.0e10])
process_nuclide("Ni61", "n-028_Ni_061.endf", [1.0e10])
process_nuclide("Ni62", "n-028_Ni_062.endf", [1.0e10])
process_nuclide("Ni64", "n-028_Ni_064.endf", [1.0e10])
process_nuclide("Cu63", "n-029_Cu_063.endf", [1.0e10])
process_nuclide("Cu65", "n-029_Cu_065.endf", [1.0e10])
process_nuclide("Br81", "n-035_Br_081.endf", [1.0e10])
process_nuclide("Kr82", "n-036_Kr_082.endf", [1.0e10])
process_nuclide("Kr83", "n-036_Kr_083.endf", [1.0e10])
process_nuclide("Kr84", "n-036_Kr_084.endf", [1.0e10])
process_nuclide("Kr85", "n-036_Kr_085.endf", [1.0e10])
process_nuclide("Kr86", "n-036_Kr_086.endf", [1.0e10])
process_nuclide("Sr89", "n-038_Sr_089.endf", [1.0e10])
process_nuclide("Sr90", "n-038_Sr_090.endf", [1.0e10])
process_nuclide("Y89", "n-039_Y_089.endf", [1.0e10])
process_nuclide("Y90", "n-039_Y_090.endf", [1.0e10])
process_nuclide("Y91", "n-039_Y_091.endf", [1.0e10])
process_nuclide("Zr90", "n-040_Zr_090.endf")
process_nuclide("Zr91", "n-040_Zr_091.endf")
process_nuclide("Zr92", "n-040_Zr_092.endf")
process_nuclide("Zr93", "n-040_Zr_093.endf", [1.0e10])  # Fission Product
process_nuclide("Zr94", "n-040_Zr_094.endf")
process_nuclide("Zr95", "n-040_Zr_095.endf", [1.0e10])  # Fission Product
process_nuclide("Zr96", "n-040_Zr_096.endf")
process_nuclide("Nb95", "n-041_Nb_095.endf", [1.0e10])  # Fission Product
process_nuclide("Mo92", "n-042_Mo_092.endf")
process_nuclide("Mo94", "n-042_Mo_094.endf")
process_nuclide("Mo95", "n-042_Mo_095.endf")
process_nuclide("Mo96", "n-042_Mo_096.endf")
process_nuclide("Mo97", "n-042_Mo_097.endf")
process_nuclide("Mo98", "n-042_Mo_098.endf")
process_nuclide("Mo99", "n-042_Mo_099.endf", [1.0e10])  # Fission Product
process_nuclide("Mo100", "n-042_Mo_100.endf")
process_nuclide("Tc99", "n-043_Tc_099.endf")
process_nuclide("Ru99", "n-044_Ru_099.endf", [1.0e10])  # Fission Product
process_nuclide("Ru100", "n-044_Ru_100.endf", [1.0e10])  # Fission Product
process_nuclide("Ru101", "n-044_Ru_101.endf", [1.0e10])  # Fission Product
process_nuclide("Ru102", "n-044_Ru_102.endf", [1.0e10])  # Fission Product
process_nuclide("Ru103", "n-044_Ru_103.endf", [1.0e10])  # Fission Product
process_nuclide("Ru104", "n-044_Ru_104.endf", [1.0e10])  # Fission Product
process_nuclide("Ru105", "n-044_Ru_105.endf", [1.0e10])  # Fission Product
process_nuclide("Ru106", "n-044_Ru_106.endf", [1.0e10])  # Fission Product
process_nuclide("Rh103", "n-045_Rh_103.endf", [1.0e10])  # Fission Product
process_nuclide("Rh105", "n-045_Rh_105.endf", [1.0e10])  # Fission Product
process_nuclide("Pd104", "n-046_Pd_104.endf", [1.0e10])  # Fission Product
process_nuclide("Pd105", "n-046_Pd_105.endf", [1.0e10])  # Fission Product
process_nuclide("Pd106", "n-046_Pd_106.endf", [1.0e10])  # Fission Product
process_nuclide("Pd107", "n-046_Pd_107.endf", [1.0e10])  # Fission Product
process_nuclide("Pd108", "n-046_Pd_108.endf", [1.0e10])  # Fission Product
process_nuclide("Ag107", "n-047_Ag_107.endf")
process_nuclide("Ag109", "n-047_Ag_109.endf")
process_nuclide("Ag110m1", "n-047_Ag_110m1.endf")
process_nuclide("Ag111", "n-047_Ag_111.endf")
process_nuclide("Cd106", "n-048_Cd_106.endf")
process_nuclide("Cd108", "n-048_Cd_108.endf")
process_nuclide("Cd110", "n-048_Cd_110.endf")
process_nuclide("Cd111", "n-048_Cd_111.endf")
process_nuclide("Cd112", "n-048_Cd_112.endf")
process_nuclide("Cd113", "n-048_Cd_113.endf")
process_nuclide("Cd114", "n-048_Cd_114.endf")
process_nuclide("Cd116", "n-048_Cd_116.endf")
process_nuclide("In113", "n-049_In_113.endf")
process_nuclide("In115", "n-049_In_115.endf")
process_nuclide("Sn112", "n-050_Sn_112.endf", [1.0e10])
process_nuclide("Sn114", "n-050_Sn_114.endf")
process_nuclide("Sn115", "n-050_Sn_115.endf")
process_nuclide("Sn116", "n-050_Sn_116.endf")
process_nuclide("Sn117", "n-050_Sn_117.endf")
process_nuclide("Sn118", "n-050_Sn_118.endf")
process_nuclide("Sn119", "n-050_Sn_119.endf")
process_nuclide("Sn120", "n-050_Sn_120.endf")
process_nuclide("Sn122", "n-050_Sn_122.endf")
process_nuclide("Sn124", "n-050_Sn_124.endf")
process_nuclide("Sb121", "n-051_Sb_121.endf", [1.0e10])  # Fission Product
process_nuclide("Sb123", "n-051_Sb_123.endf", [1.0e10])  # Fission Product
process_nuclide("Sb125", "n-051_Sb_125.endf", [1.0e10])  # Fission Product
process_nuclide("Te127m1", "n-052_Te_127m1.endf", [1.0e10])  # Fission Product
process_nuclide("Te129m1", "n-052_Te_129m1.endf", [1.0e10])  # Fission Product
process_nuclide("Te132", "n-052_Te_132.endf", [1.0e10])  # Fission Product
process_nuclide("I127", "n-053_I_127.endf")
process_nuclide("I129", "n-053_I_129.endf")
process_nuclide("I131", "n-053_I_131.endf", [1.0e10])  # Fission Product
process_nuclide("I135", "n-053_I_135.endf", [1.0e10])  # Fission Product
process_nuclide("Xe128", "n-054_Xe_128.endf", [1.0e10])  # Fission Product
process_nuclide("Xe129", "n-054_Xe_129.endf")
process_nuclide("Xe130", "n-054_Xe_130.endf", [1.0e10])  # Fission Product
process_nuclide("Xe131", "n-054_Xe_131.endf")
process_nuclide("Xe132", "n-054_Xe_132.endf", [1.0e10])  # Fission Product
process_nuclide("Xe133", "n-054_Xe_133.endf", [1.0e10])  # Fission Product
process_nuclide("Xe134", "n-054_Xe_134.endf", [1.0e10])  # Fission Product
process_nuclide("Xe135", "n-054_Xe_135.endf", [1.0e10])  # Fission Product
process_nuclide("Xe136", "n-054_Xe_136.endf", [1.0e10])  # Fission Product
process_nuclide("Cs133", "n-055_Cs_133.endf", [1.0e10])  # Fission Product
process_nuclide("Cs134", "n-055_Cs_134.endf", [1.0e10])  # Fission Product
process_nuclide("Cs135", "n-055_Cs_135.endf", [1.0e10])  # Fission Product
process_nuclide("Cs136", "n-055_Cs_136.endf", [1.0e10])  # Fission Product
process_nuclide("Cs137", "n-055_Cs_137.endf", [1.0e10])  # Fission Product
process_nuclide("Ba134", "n-056_Ba_134.endf", [1.0e10])  # Fission Product
process_nuclide("Ba137", "n-056_Ba_137.endf", [1.0e10])  # Fission Product
process_nuclide("Ba140", "n-056_Ba_140.endf", [1.0e10])  # Fission Product
process_nuclide("La139", "n-057_La_139.endf", [1.0e10])  # Fission Product
process_nuclide("La140", "n-057_La_140.endf", [1.0e10])  # Fission Product
process_nuclide("Ce140", "n-058_Ce_140.endf", [1.0e10])  # Fission Product
process_nuclide("Ce141", "n-058_Ce_141.endf", [1.0e10])  # Fission Product
process_nuclide("Ce142", "n-058_Ce_142.endf", [1.0e10])  # Fission Product
process_nuclide("Ce143", "n-058_Ce_143.endf", [1.0e10])  # Fission Product
process_nuclide("Ce144", "n-058_Ce_144.endf", [1.0e10])  # Fission Product
process_nuclide("Pr141", "n-059_Pr_141.endf", [1.0e10])  # Fission Product
process_nuclide("Pr143", "n-059_Pr_143.endf", [1.0e10])  # Fission Product
process_nuclide("Nd142", "n-060_Nd_142.endf", [1.0e10])  # Fission Product
process_nuclide("Nd143", "n-060_Nd_143.endf", [1.0e10])  # Fission Product
process_nuclide("Nd144", "n-060_Nd_144.endf", [1.0e10])  # Fission Product
process_nuclide("Nd145", "n-060_Nd_145.endf", [1.0e10])  # Fission Product
process_nuclide("Nd146", "n-060_Nd_146.endf", [1.0e10])  # Fission Product
process_nuclide("Nd147", "n-060_Nd_147.endf", [1.0e10])  # Fission Product
process_nuclide("Nd148", "n-060_Nd_148.endf", [1.0e10])  # Fission Product
process_nuclide("Nd150", "n-060_Nd_150.endf", [1.0e10])  # Fission Product
process_nuclide("Pm147", "n-061_Pm_147.endf", [1.0e10])  # Fission Product
process_nuclide("Pm148", "n-061_Pm_148.endf", [1.0e10])  # Fission Product
process_nuclide("Pm148m1", "n-061_Pm_148m1.endf", [1.0e10])  # Fission Product
process_nuclide("Pm149", "n-061_Pm_149.endf", [1.0e10])  # Fission Product
process_nuclide("Pm151", "n-061_Pm_151.endf", [1.0e10])  # Fission Product
process_nuclide("Sm147", "n-062_Sm_147.endf")
process_nuclide("Sm148", "n-062_Sm_148.endf", [1.0e10])  # Fission Product
process_nuclide("Sm149", "n-062_Sm_149.endf")
process_nuclide("Sm150", "n-062_Sm_150.endf")
process_nuclide("Sm151", "n-062_Sm_151.endf")
process_nuclide("Sm152", "n-062_Sm_152.endf")
process_nuclide("Sm153", "n-062_Sm_153.endf")
process_nuclide("Sm154", "n-062_Sm_154.endf")
process_nuclide("Eu151", "n-063_Eu_151.endf")
process_nuclide("Eu152", "n-063_Eu_152.endf")
process_nuclide("Eu153", "n-063_Eu_153.endf")
process_nuclide("Eu154", "n-063_Eu_154.endf")
process_nuclide("Eu155", "n-063_Eu_155.endf")
process_nuclide("Eu156", "n-063_Eu_156.endf", [1.0e10])  # Fission Product
process_nuclide("Gd152", "n-064_Gd_152.endf")
process_nuclide("Gd154", "n-064_Gd_154.endf")
process_nuclide("Gd155", "n-064_Gd_155.endf")
process_nuclide("Gd156", "n-064_Gd_156.endf")
process_nuclide("Gd157", "n-064_Gd_157.endf")
process_nuclide("Gd158", "n-064_Gd_158.endf")
process_nuclide("Gd160", "n-064_Gd_160.endf")
process_nuclide("Tb159", "n-065_Tb_159.endf", [1.0e10])  # Fission Product
process_nuclide("Tb160", "n-065_Tb_160.endf", [1.0e10])  # Fission Product
process_nuclide("Tb161", "n-065_Tb_161.endf", [1.0e10])  # Fission Product
process_nuclide("Dy160", "n-066_Dy_160.endf")
process_nuclide("Dy161", "n-066_Dy_161.endf")
process_nuclide("Dy162", "n-066_Dy_162.endf")
process_nuclide("Dy163", "n-066_Dy_163.endf")
process_nuclide("Dy164", "n-066_Dy_164.endf")
process_nuclide("Ho165", "n-067_Ho_165.endf", [1.0e10])  # Fission Product
process_nuclide("Er162", "n-068_Er_162.endf", [1.0e10])  # Fission Product
process_nuclide("Er164", "n-068_Er_164.endf", [1.0e10])  # Fission Product
process_nuclide("Er166", "n-068_Er_166.endf", [1.0e10])  # Fission Product
process_nuclide("Er167", "n-068_Er_167.endf", [1.0e10])  # Fission Product
process_nuclide("Er168", "n-068_Er_168.endf", [1.0e10])  # Fission Product
process_nuclide("Er169", "n-068_Er_169.endf", [1.0e10])  # Fission Product
process_nuclide("Er170", "n-068_Er_170.endf", [1.0e10])  # Fission Product
process_nuclide("Tm169", "n-069_Tm_169.endf", [1.0e10])  # Fission Product
process_nuclide("Tm170", "n-069_Tm_170.endf", [1.0e10])  # Fission Product
process_nuclide("Tm171", "n-069_Tm_171.endf", [1.0e10])  # Fission Product
process_nuclide("Hf174", "n-072_Hf_174.endf")
process_nuclide("Hf176", "n-072_Hf_176.endf")
process_nuclide("Hf177", "n-072_Hf_177.endf")
process_nuclide("Hf178", "n-072_Hf_178.endf")
process_nuclide("Hf179", "n-072_Hf_179.endf")
process_nuclide("Hf180", "n-072_Hf_180.endf")
process_nuclide("Hf181", "n-072_Hf_181.endf")
process_nuclide("Ta181", "n-073_Ta_181.endf")
process_nuclide("Ta182", "n-073_Ta_182.endf", [1.0e10])  # Fission Product
process_nuclide("Th230", "n-090_Th_230.endf")
process_nuclide("Th231", "n-090_Th_231.endf")
process_nuclide("Th232", "n-090_Th_232.endf")
process_nuclide("Th234", "n-090_Th_234.endf")
process_nuclide("Pa231", "n-091_Pa_231.endf")
process_nuclide("Pa232", "n-091_Pa_232.endf")
process_nuclide("Pa233", "n-091_Pa_233.endf")
process_nuclide("U232", "n-092_U_232.endf")
process_nuclide("U233", "n-092_U_233.endf")
process_nuclide("U234", "n-092_U_234.endf")
process_nuclide("U236", "n-092_U_236.endf")
process_nuclide("U237", "n-092_U_237.endf")
process_nuclide("U238", "n-092_U_238.endf")
process_nuclide("Np236", "n-093_Np_236.endf")
process_nuclide("Np237", "n-093_Np_237.endf")
process_nuclide("Np238", "n-093_Np_238.endf")
process_nuclide("Np239", "n-093_Np_239.endf")
process_nuclide("Pu236", "n-094_Pu_236.endf")
process_nuclide("Pu237", "n-094_Pu_237.endf")
process_nuclide("Pu238", "n-094_Pu_238.endf")
process_nuclide("Pu239", "n-094_Pu_239.endf")
process_nuclide("Pu240", "n-094_Pu_240.endf")
process_nuclide("Pu241", "n-094_Pu_241.endf")
process_nuclide("Pu242", "n-094_Pu_242.endf")
process_nuclide("Pu243", "n-094_Pu_243.endf")
process_nuclide("Pu244", "n-094_Pu_244.endf")
process_nuclide("Am241", "n-095_Am_241.endf")
process_nuclide("Am242m1", "n-095_Am_242m1.endf")
process_nuclide("Am243", "n-095_Am_243.endf")
process_nuclide("Cm242", "n-096_Cm_242.endf")
process_nuclide("Cm243", "n-096_Cm_243.endf")
process_nuclide("Cm244", "n-096_Cm_244.endf")
process_nuclide("Cm245", "n-096_Cm_245.endf")
process_nuclide("Cm246", "n-096_Cm_246.endf")

h5.close()
