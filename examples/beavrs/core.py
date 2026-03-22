from scarabee.coeur import SimpleTile, QuadrantsTile, CoreBuilder
import matplotlib.pyplot as plt
import numpy as np
import pickle

a1_00_ = pickle.loads(open("F16_0.pkl", "rb").read())
a2_00_ = pickle.loads(open("F24_0.pkl", "rb").read())
a3_00_ = pickle.loads(open("F31_0.pkl", "rb").read())
a2_12_ = pickle.loads(open("F24_12.pkl", "rb").read())
a2_16_ = pickle.loads(open("F24_16.pkl", "rb").read())

a3_06U  = pickle.loads(open("F31_6R.pkl", "rb").read())
a3_06U.rotate_counterclockwise()

a3_06R  = pickle.loads(open("F31_6R.pkl", "rb").read())

a3_06D  = pickle.loads(open("F31_6R.pkl", "rb").read())
a3_06D.rotate_clockwise()

a3_06L  = pickle.loads(open("F31_6R.pkl", "rb").read())
a3_06L.rotate_clockwise().rotate_clockwise()

a3_152  = pickle.loads(open("F31_15II.pkl", "rb").read())

a3_151  = pickle.loads(open("F31_15II.pkl", "rb").read())
a3_151.rotate_clockwise()

a3_154  = pickle.loads(open("F31_15II.pkl", "rb").read())
a3_154.rotate_clockwise().rotate_clockwise()

a3_153  = pickle.loads(open("F31_15II.pkl", "rb").read())
a3_153.rotate_clockwise().rotate_clockwise().rotate_clockwise()

a3_16_  = pickle.loads(open("F31_16.pkl", "rb").read())
a3_20_  = pickle.loads(open("F31_20.pkl", "rb").read())
rf____  = pickle.loads(open("reflector.pkl", "rb").read())

# Define core geometry
#                R       P       N       M       L       K       J       H       G       F       E       D       C       B       A
tiles = [[[0.    , 0.    , 0.    , 0.    , rf____, rf____, rf____, rf____, rf____, rf____, rf____, rf____, rf____, 0.    , 0.    , 0.    , 0.],

          [0.    , 0.    , rf____, rf____, rf____, a3_00_, a3_06D, a3_00_, a3_06D, a3_00_, a3_06D, a3_00_, rf____, rf____, rf____, 0.    , 0.],     #  1 

          [0.    , rf____, rf____, a3_00_, a3_00_, a3_16_, a1_00_, a3_20_, a1_00_, a3_20_, a1_00_, a3_16_, a3_00_, a3_00_, rf____, rf____, 0.],     #  2 

          [0.    , rf____, a3_00_, a3_154, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a3_153, a3_00_, rf____, 0.],     #  3

          [rf____, rf____, a3_00_, a2_16_, a2_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a2_00_, a2_16_, a3_00_, rf____, rf____], #  4

          [rf____, a3_00_, a3_16_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_16_, a3_00_, rf____], #  5

          [rf____, a3_06R, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_06L, rf____], #  6

          [rf____, a3_00_, a3_20_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a3_20_, a3_00_, rf____], #  7

          [rf____, a3_06R, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_06L, rf____], #  8

          [rf____, a3_00_, a3_20_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a3_20_, a3_00_, rf____], #  9

          [rf____, a3_06R, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_06L, rf____], # 10

          [rf____, a3_00_, a3_16_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_16_, a3_00_, rf____], # 11

          [rf____, rf____, a3_00_, a2_16_, a2_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a2_00_, a2_16_, a3_00_, rf____, rf____], # 12

          [0.    , rf____, a3_00_, a3_151, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a3_152, a3_00_, rf____, 0.],     # 13

          [0.    , rf____, rf____, a3_00_, a3_00_, a3_16_, a1_00_, a3_20_, a1_00_, a3_20_, a1_00_, a3_16_, a3_00_, a3_00_, rf____, rf____, 0.],     # 14

          [0.    , 0.    , rf____, rf____, rf____, a3_00_, a3_06U, a3_00_, a3_06U, a3_00_, a3_06U, a3_00_, rf____, rf____, rf____, 0.    , 0.],     # 15

          [0.    , 0.    , 0.    , 0.    , rf____, rf____, rf____, rf____, rf____, rf____, rf____, rf____, rf____, 0.    , 0.    , 0.    , 0.],
        ]]
tiles = np.array(tiles)
dz = np.array([20.])

core_builder = CoreBuilder(21.50364, 17, 1.25984, 17, tiles, dz, 1., 1.)
core_builder.solver.solve()

# =============================================================================
# Rasterize and plot the flux and homogeneous power

# x, y, z arrays for rasterizing the flux / power
x_max = np.sum(core_builder.dx)
dx = x_max / 500.
x = np.arange(start=0., stop=x_max+dx, step=dx)
y_max = np.sum(core_builder.dy)
dy = y_max / 500.
y = np.arange(start=0., stop=y_max+dy, step=dy)
z = np.array([0.])

# Plot the flux in each group
flux = core_builder.solver.flux(x, y, z)[:,:,:,0]
for g in range(flux.shape[0]):
    plt.pcolormesh(y, x, flux[g,:,:], cmap='turbo')
    plt.title(f"Group {g+1} Flux")
    plt.show()

# Plot the homogeneous power
power = core_builder.solver.power(x, y, z)[:,:,0]
plt.pcolormesh(y, x, power, cmap='turbo')
plt.title("Homogeneous Power Distribution")
plt.show()

# =============================================================================
# Assembly Powers

asmbly_powers = core_builder.compute_assembly_powers()
plt.imshow(asmbly_powers, cmap='turbo')
plt.title("Assembly Powers")
plt.show()

# =============================================================================
# Pin Powers
pin_power = core_builder.compute_pin_powers(np.array([10.]))
print("Max Pin Power: {:.5f}".format(np.nanmax(pin_power)))
print("Min Pin Power: {:.5f}".format(np.nanmin(pin_power[np.where(pin_power != 0.)])))
pin_power[np.where(pin_power == 0.)] = np.nan

plt.pcolormesh(core_builder.x_pin_centers, core_builder.y_pin_centers, pin_power[:,:,0], cmap='turbo')
plt.title("Pin Power Distribution")
plt.show()
