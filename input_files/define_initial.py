itype = 5

if itype == 1:
  # Guassian in frequency and space, cos type in direction
  fp = 0.10
  fp_spread = 0.01
  mean_direction = 270
  cosine_power = 2
  Xm = 0.0
  X_spread = 1000.0
  Ym = 1.0
  Y_spread = 1000.0
  Hmax = 2.5

elif itype == 2:
  # JONSWAP with Hassselman et al. (1980) direct. distribution
  alpha = 0.0081
  peak_freq = 0.1
  mean_direction = 270
  gamma = 1.0
  sigA = 0.0
  sigB = 0.0
  Xm = 1.0
  X_spread = 100.0
  Ym = 1.0
  Y_spread = 100.0

elif itype == 3:
  # Fetch limited JONSWAP
  pass
elif itype == 4:
  # User-defined spectrum
  scale_factor = -0.1
  F = []
elif itype == 5:
  # Start from calm conditions
  pass
