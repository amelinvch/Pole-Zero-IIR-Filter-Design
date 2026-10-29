from pathlib import Path
import numpy as np
from scipy.io.wavfile import write

fs = 4000
duration = 10
t = np.arange(fs * duration) / fs
rng = np.random.default_rng(42)
output = Path(__file__).resolve().parent / "signals"
output.mkdir(exist_ok=True)

# Synthetic ECG-like waveform (not a medical recording)
ecg_clean = np.zeros(len(t))
for beat in np.arange(0.5, duration, 1.0):
  ecg_clean += 0.12 * np.exp(-0.5 * ((t - (beat - 0.20)) / 0.045)**2)
  ecg_clean -= 0.15 * np.exp(-0.5 * ((t - (beat - 0.025)) / 0.010)**2)
  ecg_clean += 1.00 * np.exp(-0.5 * ((t - beat) / 0.012)**2)
  ecg_clean -= 0.25 * np.exp(-0.5 * ((t - (beat + 0.025)) / 0.012)**2)
  ecg_clean += 0.30 * np.exp(-0.5 * ((t - (beat + 0.28)) / 0.080)**2)
ecg_noisy = ecg_clean + 0.04 * np.sin(2 * np.pi * 0.3 * t)
ecg_noisy += 0.20 * np.sin(2 * np.pi * 50 * t)
ecg_noisy += 0.015 * rng.standard_normal(len(t))

# Synthetic mixture; clean reference is the desired 100 Hz component
envelope = 0.5 + 0.5 * np.sin(2 * np.pi * 0.4 * t)**2
tone_clean = 0.08 * envelope * np.sin(2 * np.pi * 100 * t)
tone_noisy = tone_clean + 0.20 * np.sin(2 * np.pi * 40 * t)
tone_noisy += 0.15 * np.sin(2 * np.pi * 180 * t)
tone_noisy += 0.10 * np.sin(2 * np.pi * 300 * t)
tone_noisy += 0.04 * rng.standard_normal(len(t))

for name, x in [("ecg_clean", ecg_clean), ("ecg_noisy", ecg_noisy), ("tone_clean", tone_clean), ("tone_noisy", tone_noisy)]:
  write(output / f"{name}.wav", fs, x.astype(np.float32))
  print("Created", name + ".wav")
