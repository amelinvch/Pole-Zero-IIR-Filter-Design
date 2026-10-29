import math

def design_biquad_filter(G0, Gc, B, wc):
  if not (Gc > 0 and G0 > 0):
    raise ValueError("G0 and Gc must be positive")
  if not (0 < wc < math.pi):
    raise ValueError("Center frequency must be in the range (0, pi)")
  if not (0 < B < math.pi):
    raise ValueError("Bandwidth must be in the range (0, pi)")

  gain_ratio = Gc / G0

  if gain_ratio == 1:
    return [G0, 0, 0], [1, 0, 0]
  
  if not (0 < gain_ratio < 1 or gain_ratio >= 2):
    raise ValueError("gain_ratio must be in the range (0, 1) or [2, +inf)")

  if (0 < gain_ratio < 1): # notch_filter
    if (0 < gain_ratio <= 0.5):
      rz = 1 - (B / 2) * math.sqrt(1 - 2 * gain_ratio**2)
      rp = 1 - (B / 2) * math.sqrt(1/gain_ratio**2 - 2)
    elif (0.5 < gain_ratio < 1):
      rz = 1 - (B / 2) * math.sqrt(gain_ratio)
      rp = 1 - B / (2 * math.sqrt(gain_ratio))
    K = G0

  elif (gain_ratio >= 2): # resonance_filter
    rp = 1 - B / 2
    val = (B * gain_ratio / 2) * (math.sin(wc) / (1 - math.cos(wc)))
    rz = (1 - val) / ( 1- val + B * gain_ratio / 2)
    K = G0 * (1 - 2 * rp * math.cos(wc) + rp**2) / (1 - 2 * rz * math.cos(wc) + rz**2)

  b = [K, -2 * K * rz * math.cos(wc), K * rz**2]
  a = [1, -2 * rp * math.cos(wc), rp**2]

  return b, a