"""
Joint design parameters and calibrated model coefficients.
Source: Gaozhang et al. (2024), Table 1 and Section 2.
All values verified against the paper.
"""

# Joint geometry
JOINT_HEIGHT_M = 0.049        # 49 mm
JOINT_DIAMETER_M = 0.038      # 38 mm
CHAMBER_SHAPE = "semi-circular"

# Design parameters (Table 1)
n0 = 4.601e-10
Ap0 = 2.395e-4    # m^2, initial cavity cross-section area
Ai = 1.793e-4     # m^2, chamber wall cross-section area
a = 0.018          # m, outer radius
u_bar = 0.207e6    # Pa (0.207 MPa), effective initial shear modulus
n1 = 100.844
n2 = 1.4
Lm = 0.009         # m, moment arm
L = 0.049          # m, chamber length

# Stiffness model Fourier coefficients (Table 1)
a0_coeff = 7.542
a1_coeff = -2.183
a2_coeff = -1.304
b0_coeff = 0.054
b1_coeff = -4.935
b2_coeff = 0.427
m1_coeff = 1.159e-4
m2_coeff = -0.129
w_coeff = 0.067

# Operating range
MAX_PRESSURE_PA = 3.0e5       # 3 x 10^5 Pa
PRESSURE_STEP_PA = 0.5e5      # 0.5 x 10^5 Pa

# Key experimental results (from paper text and figures)
MAX_BENDING_ANGLE_DEG = 48.8  # one direction, at P1=3e5, P2=0
TOTAL_BENDING_RANGE_DEG = 100 # approximate, both directions
MAX_OUTPUT_FORCE_N = 20.0     # at 1.5e5 Pa
STIFFNESS_MIN = 26.56         # N.mm/deg, at P2=0, P1=1e5
STIFFNESS_MAX = 102.59        # N.mm/deg, at P1=P2=3e5

# Paper's own model accuracy
KINEMATIC_MODEL_AVG_DEVIATION_DEG = 0.927
KINEMATIC_MODEL_AVG_DEVIATION_RATE = 0.0193  # 1.93%
STIFFNESS_MODEL_AVG_DEVIATION = 5.394        # N.mm/deg
STIFFNESS_MODEL_DEVIATION_RATE = 0.0525      # 5.25%
