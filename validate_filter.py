
import math
import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import freqz
from design_biquad_filter import design_biquad_filter


# General parameters
fs = 4000

# (name, center frequency, bandwidth, G0, Gc)
cases = [
    ("Notch - 5 Hz", 50, 5, 1, 0.1),
    ("Notch - 0.5 Hz", 50, 0.5, 1, 0.1),
    ("Resonance - 10 Hz", 100, 10, 1, 5),
    ("Resonance - 1 Hz", 100, 1, 1, 5)
]


fig, axes = plt.subplots(2, 2, figsize=(13, 9))
axes = axes.flatten()


for i in range(len(cases)):
    name, fc, bandwidth, G0, Gc = cases[i]

    # Convert Hz to rad/sample
    wc = 2 * math.pi * fc / fs
    B = 2 * math.pi * bandwidth / fs

    # Design filter
    b, a = design_biquad_filter(G0, Gc, B, wc)

    # Frequency response
    frequencies, h = freqz(b, a, worN=131072, fs=fs)
    magnitude = np.abs(h)

    # Calculate gains
    dc_gain = abs(sum(b) / sum(a))

    _, h_center = freqz(b, a, worN=[wc])
    center_gain = abs(h_center[0])

    # Find actual center frequency
    if Gc < G0:
        index = np.argmin(magnitude)
    else:
        index = np.argmax(magnitude)

    actual_fc = frequencies[index]
    extremum_gain = magnitude[index]

    # Calculate errors
    dc_error = 100 * (dc_gain - G0) / G0
    gain_error = 100 * (center_gain - Gc) / Gc
    frequency_error = actual_fc - fc

    # Check stability
    poles = np.roots(a)
    zeros = np.roots(b)

    pole_radius = abs(poles[0])
    zero_radius = abs(zeros[0])

    stable = np.all(np.abs(poles) < 1)

    if Gc < G0:
        placement_valid = pole_radius < zero_radius <= 1
    else:
        placement_valid = 0 <= zero_radius < pole_radius < 1

    # Bandwidth level (Moorer's definition)
    gain_ratio = Gc / G0

    if gain_ratio <= 0.5:
        level = extremum_gain * math.sqrt(2)
    elif gain_ratio < 1:
        level = math.sqrt(extremum_gain * dc_gain)
    else:
        level = extremum_gain / math.sqrt(2)

    # Find bandwidth edges
    lower_frequency = None
    upper_frequency = None

    for j in range(len(frequencies) - 1):

        y1 = magnitude[j]
        y2 = magnitude[j + 1]

        # Check if the response crosses the bandwidth level
        if (y1 - level) * (y2 - level) < 0:

            f1 = frequencies[j]
            f2 = frequencies[j + 1]

            # Linear interpolation
            crossing = f1 + (level - y1) * (f2 - f1) / (y2 - y1)

            if j < index:
                lower_frequency = crossing

            elif upper_frequency is None:
                upper_frequency = crossing

    # Calculate bandwidth
    actual_bandwidth = None
    bandwidth_error = None

    if lower_frequency is not None and upper_frequency is not None:
        actual_bandwidth = upper_frequency - lower_frequency
        bandwidth_error = 100 * (actual_bandwidth - bandwidth) / bandwidth

    # Print results
    print("\n" + "=" * 45)
    print(name)
    print("=" * 45)

    print(f"Target center frequency: {fc} Hz")
    print(f"Target bandwidth: {bandwidth} Hz")
    print(f"G0 = {G0}, Gc = {Gc}")

    print("\nCoefficients:")
    print("b =", b)
    print("a =", a)

    print("\nValidation:")
    print(f"DC gain: {dc_gain:.6f} (error: {dc_error:+.3f}%)")
    print(f"Center gain: {center_gain:.6f} (error: {gain_error:+.3f}%)")
    print(f"Center frequency: {actual_fc:.6f} Hz (error: {frequency_error:+.6f} Hz)")

    if actual_bandwidth is not None:
        print(f"Bandwidth: {actual_bandwidth:.6f} Hz (error: {bandwidth_error:+.3f}%)")
        print(f"Lower edge: {lower_frequency:.6f} Hz")
        print(f"Upper edge: {upper_frequency:.6f} Hz")
    else:
        print("Bandwidth could not be determined")

    print("Stable:", stable)
    print("Pole-zero placement valid:", placement_valid)
    print("Pole radius:", pole_radius)
    print("Zero radius:", zero_radius)

    # Plot frequency response
    ax = axes[i]

    ax.plot(frequencies, magnitude, label="Frequency response")
    ax.axvline(fc, color="red", linestyle="--", label="Target frequency")
    ax.axhline(Gc, color="green", linestyle=":", label="Target gain")

    if actual_bandwidth is not None:
        ax.axvline(lower_frequency, color="orange", linestyle="--", label="Bandwidth edges")
        ax.axvline(upper_frequency, color="orange", linestyle="--")
        ax.axhline(level, color="purple", linestyle=":", label="Bandwidth level")

    ax.set_xlim(
        max(0, fc - 5 * bandwidth),
        min(fs / 2, fc + 5 * bandwidth)
    )

    ax.set_title(name)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude")
    ax.grid(True)
    ax.legend()


plt.tight_layout()
plt.show()
