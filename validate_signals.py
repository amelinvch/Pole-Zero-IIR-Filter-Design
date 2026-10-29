from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.io.wavfile import read
from scipy.signal import freqz, lfilter

from design_biquad_filter import design_biquad_filter


folder = Path(__file__).resolve().parent
signal_folder = folder / "signals"
figure_folder = folder / "images"
figure_folder.mkdir(exist_ok=True)

# All signals were generated at the same sampling rate.
fs = 4000

# Signal 1: ECG-like waveform with 50 Hz power-line interference.
fs_ecg, ecg_clean = read(signal_folder / "ecg_clean.wav")
fs_noisy_ecg, ecg_noisy = read(signal_folder / "ecg_noisy.wav")

# Signal 2: amplitude-modulated 100 Hz tone with other components.
fs_tone, tone_clean = read(signal_folder / "tone_clean.wav")
fs_noisy_tone, tone_noisy = read(signal_folder / "tone_noisy.wav")

if not all(rate == fs for rate in [fs_ecg, fs_noisy_ecg, fs_tone, fs_noisy_tone]):
    raise ValueError("All WAV files must have a sampling rate of 4000 Hz")

# Frequency parameters in Hz are converted to radians per sample.
notch_fc = 50
notch_bw = 5
notch_G0 = 1
notch_Gc = 0.1

resonance_fc = 100
resonance_bw = 10
resonance_G0 = 1
resonance_Gc = 5

notch_wc = 2 * math.pi * notch_fc / fs
notch_B = 2 * math.pi * notch_bw / fs
resonance_wc = 2 * math.pi * resonance_fc / fs
resonance_B = 2 * math.pi * resonance_bw / fs

b_notch, a_notch = design_biquad_filter(notch_G0, notch_Gc, notch_B, notch_wc)
b_res, a_res = design_biquad_filter(
    resonance_G0, resonance_Gc, resonance_B, resonance_wc
)

# Filter the noisy signals. Use the same implementation for both cases.
ecg_filtered = lfilter(b_notch, a_notch, ecg_noisy)
tone_filtered = lfilter(b_res, a_res, tone_noisy)

for name, b, a, wc, target in [
    ("Notch", b_notch, a_notch, notch_wc, notch_Gc),
    ("Resonance", b_res, a_res, resonance_wc, resonance_Gc),
]:
    poles = np.roots(a)
    _, response = freqz(b, a, worN=[wc])
    actual_gain = abs(response[0])
    print(f"{name}: stable = {np.all(np.abs(poles) < 1)}")
    print(f"{name}: target center gain = {target:.4f}, "
          f"actual center gain = {actual_gain:.4f}")
    if not np.all(np.abs(poles) < 1):
        raise ValueError(f"{name} filter is unstable")

# ECG: the clean reference allows a quantitative comparison.
mse_before = np.mean((ecg_noisy - ecg_clean)**2)
mse_after = np.mean((ecg_filtered - ecg_clean)**2)
mse_reduction = 100 * (mse_before - mse_after) / mse_before
print(f"ECG MSE before: {mse_before:.6f}")
print(f"ECG MSE after:  {mse_after:.6f}")
print(f"ECG MSE reduction: {mse_reduction:.2f}%")

# Plot the first two seconds of the ECG-like signal.
t_ecg = np.arange(len(ecg_clean)) / fs
fig, axes = plt.subplots(2, 1, figsize=(11, 7))
view_ecg = t_ecg < 2
axes[0].plot(t_ecg[view_ecg], ecg_clean[view_ecg], label="Clean", linewidth=1.5)
axes[0].plot(t_ecg[view_ecg], ecg_noisy[view_ecg], label="Noisy", alpha=0.6)
axes[0].plot(t_ecg[view_ecg], ecg_filtered[view_ecg], label="Filtered", linewidth=1)
axes[0].set(title="ECG-like signal: 50 Hz notch", xlabel="Time (s)", ylabel="Amplitude")
axes[0].grid(True)
axes[0].legend()

freq_ecg = np.fft.rfftfreq(len(ecg_noisy), 1 / fs)
fft_before = np.abs(np.fft.rfft(ecg_noisy)) / len(ecg_noisy)
fft_after = np.abs(np.fft.rfft(ecg_filtered)) / len(ecg_filtered)
axes[1].plot(freq_ecg, fft_before, label="Before filtering")
axes[1].plot(freq_ecg, fft_after, label="After filtering")
axes[1].axvline(50, color="red", linestyle="--", label="50 Hz interference")
axes[1].set(xlim=(0, 120), xlabel="Frequency (Hz)", ylabel="FFT magnitude")
axes[1].grid(True)
axes[1].legend()
fig.tight_layout()
fig.savefig(figure_folder / "ecg_notch.png", dpi=160)

# Resonance is a frequency-selective boost, not a denoising operation.
# Compare the filtered output against the input and a gain-adjusted 100 Hz reference.
t_tone = np.arange(len(tone_clean)) / fs
fig, axes = plt.subplots(3, 1, figsize=(11, 10))
view_tone = (t_tone >= 1.0) & (t_tone < 1.15)  # skip the initial filter transient
axes[0].plot(t_tone[view_tone], tone_noisy[view_tone], label="Input mixture", alpha=0.8)
axes[0].plot(t_tone[view_tone], tone_filtered[view_tone], label="Resonance output", alpha=0.8)
axes[0].set(title="100 Hz resonance: input and output (steady state)",
            xlabel="Time (s)", ylabel="Amplitude")
axes[0].grid(True)
axes[0].legend()

freq_tone = np.fft.rfftfreq(len(tone_noisy), 1 / fs)
tone_before = np.abs(np.fft.rfft(tone_noisy)) / len(tone_noisy)
tone_after = np.abs(np.fft.rfft(tone_filtered)) / len(tone_filtered)
axes[1].plot(freq_tone, tone_before, label="Input spectrum", linewidth=1.5)
axes[1].plot(freq_tone, tone_after, label="Output spectrum", linewidth=1.5)
axes[1].axvline(100, color="black", linestyle="--", alpha=0.6, label="100 Hz center")
axes[1].set(xlim=(0, 350), xlabel="Frequency (Hz)", ylabel="FFT magnitude",
            title="Full spectrum: other tones remain largely unchanged")
axes[1].grid(True)
axes[1].legend()

axes[2].plot(freq_tone, tone_before, label="Input spectrum", linewidth=1.5)
axes[2].plot(freq_tone, tone_after, label="Output spectrum", linewidth=1.5)
axes[2].axvline(100, color="black", linestyle="--", alpha=0.6)
axes[2].set(xlim=(90, 110), xlabel="Frequency (Hz)", ylabel="FFT magnitude",
            title="Zoom around 100 Hz: the resonant component is boosted")
axes[2].grid(True)
axes[2].legend()

# A narrowband energy ratio is more informative than reading a single FFT bin
# because the target tone is amplitude-modulated.
band = (freq_tone >= 98) & (freq_tone <= 102)
band_rms_before = np.sqrt(np.sum(tone_before[band] ** 2))
band_rms_after = np.sqrt(np.sum(tone_after[band] ** 2))
print(f"100 Hz band (98-102 Hz) magnitude ratio: "
      f"{band_rms_after / band_rms_before:.3f}x")
fig.tight_layout()
fig.savefig(figure_folder / "tone_resonance.png", dpi=160)

# Compare the actual filter frequency responses with the target frequencies.
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for ax, b, a, fc, Gc, title, limits in [
    (axes[0], b_notch, a_notch, notch_fc, notch_Gc, "Notch response", (30, 70)),
    (axes[1], b_res, a_res, resonance_fc, resonance_Gc, "Resonance response", (75, 125)),
]:
    freq, response = freqz(b, a, worN=131072, fs=fs)
    ax.plot(freq, np.abs(response), label="Actual magnitude")
    ax.axvline(fc, color="red", linestyle="--", label="Target frequency")
    ax.axhline(Gc, color="green", linestyle=":", label="Target gain")
    ax.set(xlim=limits, title=title, xlabel="Frequency (Hz)", ylabel="Magnitude")
    ax.grid(True)
    ax.legend()
fig.tight_layout()
fig.savefig(figure_folder / "resonance_notch_frequency_responses.png", dpi=160)

print("\nInterpretation:")
print("- Notch: the 50 Hz spectral component is reduced, as intended.")
print("- Resonance: energy around 100 Hz is amplified, as intended.")
print("- Other frequencies are not completely removed by the resonance filter.")
print("- Approximate pole-zero equations do not guarantee exact target gain.")
print("- These are synthetic examples of realistic interference scenarios.")
print("\nFigures saved in:", figure_folder)
plt.show()
