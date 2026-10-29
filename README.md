# Second-Order IIR Parametric Equalizer Filter Design

**DSP assignment — implementation and evaluation of the pole–zero placement method**

This project implements the second-order IIR filter design approach described by Van Waterschoot and Moonen (2007), *A Pole-Zero Placement Technique for Designing Second-Order IIR Parametric Equalizer Filters*. It examines two filter types: **notch** (attenuation around a chosen frequency) and **resonance** (amplification around a chosen frequency).

The aim is to implement the coefficient calculations, check the resulting frequency responses against the requested specifications, and demonstrate their effect on example signals.

## 1. Method

The filter is represented by

$$
H(z)=K\frac{1-2r_z\cos(\omega_c)z^{-1}+r_z^2z^{-2}}
{1-2r_p\cos(\omega_c)z^{-1}+r_p^2z^{-2}},
$$

where $r_z$ and $r_p$ are the zero and pole radii, $\omega_c$ is the specified center angular frequency, and $K$ is the gain factor. The coefficients returned by `design_biquad_filter()` are

$$
\begin{aligned}
b &= [K,\;-2Kr_z\cos(\omega_c),\;Kr_z^2],\\
a &= [1,\;-2r_p\cos(\omega_c),\;r_p^2].
\end{aligned}
$$

The implementation takes `G0` (reference/DC gain), `Gc` (desired center gain), `B` (bandwidth in rad/sample), and `wc` (center frequency in rad/sample). Frequency specifications in Hz are converted using $\omega_c=2\pi f_c/f_s$ and $B=2\pi\mathrm{BW}/f_s$.

The implementation uses the paper's **approximate** pole–zero placement equations. It supports attenuation cases $0<G_c/G_0<1$, amplification cases $G_c/G_0\ge 2$, and the trivial constant-gain case $G_c/G_0=1$. The intermediate amplification range $1<G_c/G_0<2$ is not implemented.

## 2. Files

| File | Purpose |
|---|---|
| `design_biquad_filter.py` | Calculates second-order IIR numerator and denominator coefficients. |
| `validate_filter.py` | Checks gain, center frequency, bandwidth, pole/zero placement, and stability for four designs; plots their responses. |
| `generate_signals.py` | Generates reproducible synthetic ECG-like and audio-like WAV signals. |
| `filter_analysis.ipynb` | Applies the filters to signals, compares waveforms and spectra, and computes evaluation metrics. |
| `Figure_1.png` | Frequency responses exported directly from `validate_filter.py`. |


## 3. Numerical validation

All four configurations use a sampling frequency of **4,000 Hz**, reference gain $G_0=1$, and the parameters below.

| Filter | Center frequency | Target bandwidth | Target center gain | Actual center gain | Center-gain error | Measured bandwidth | BW error |
|---|---:|---:|---:|---:|---:|---:|---:|
| Notch | 50 Hz | 5 Hz | 0.1 | 0.098732 | −1.27% | 5.008561 Hz | +0.17% |
| Notch | 50 Hz | 0.5 Hz | 0.1 | 0.100145 | +0.15% | 0.500070 Hz | +0.01% |
| Resonance | 100 Hz | 10 Hz | 5 | 8.014332 | +60.29% | 10.152023 Hz | +1.52% |
| Resonance | 100 Hz | 1 Hz | 5 | 5.246629 | +4.93% | 1.039014 Hz | +3.90% |

The response is evaluated with `scipy.signal.freqz`. The measured bandwidth is the distance between the two interpolated crossings of the bandwidth level around the response minimum/maximum. The bandwidth levels follow the definitions used in the validation script (including the $\sqrt{2}$ threshold for the tested notch and resonance cases).

**Figure 1.** Frequency responses of the four tested IIR filter configurations. The blue curves show the measured magnitude response, red dashed lines mark the target center frequencies, green dotted lines mark the desired center gains, orange dashed lines show the measured bandwidth boundaries, and purple dotted lines indicate the bandwidth measurement levels.

![Figure 1 — Four filter frequency responses](Figure_1.png)

**Stability and placement.** All four tested designs have poles strictly inside the unit circle. Their pole/zero radii also satisfy the expected ordering: poles closer to the origin than zeros for the tested notch filters, and zeros closer to the origin than poles for the tested resonance filters. These checks confirm stability for the four selected parameter sets, not for every possible input.

**Reference gain.** The measured DC gains for the two notch cases are **0.827824** (5 Hz bandwidth) and **1.001074** (0.5 Hz bandwidth). The first is **17.22% below** the specified $G_0=1$. Both tested resonance filters have DC gain approximately equal to 1.

## 4. Interpretation of the results

The **0.5 Hz notch** closely matches the target center gain and bandwidth. The **5 Hz notch** also gives a good center-gain and bandwidth match, but its DC gain is noticeably lower than requested.

The **10 Hz resonance** exhibits the largest discrepancy: the center gain is about **8.01 instead of 5**, an error of **+60.29%**. Reducing the resonance bandwidth to **1 Hz** brings the center gain closer to the target (**5.25**, +4.93%), although its relative bandwidth error is larger than in the 10 Hz case.

These results demonstrate a limitation of the **approximate design equations**: satisfying the pole/zero geometry and producing a stable filter does **not** guarantee that every specified gain and bandwidth will be matched accurately. The observed trends apply to these test configurations; they are not a general proof that narrower bandwidth always improves every metric.

The reported center-frequency deviations are only a few thousandths of a hertz and are measured on a discrete `freqz` grid; they should not be interpreted as exact continuous-frequency peak shifts.

## 5. Signal experiments

## Signal Experiments

To evaluate the practical behaviour of the implemented filters, two synthetic signals were generated using `generate_signals.py`. Both signals were sampled at 4000 Hz and had a duration of 10 seconds.

The first signal represents an ECG-like waveform contaminated with 50 Hz power-line interference, low-frequency baseline wander and white noise. The second represents an amplitude-modulated 100 Hz tone mixed with additional frequency components and white noise.

The signals are synthetic rather than real recordings. This allows comparison against known clean references, but does not establish performance on real-world measurements.

### ECG-like Signal — Notch Filtering

A notch filter was applied to suppress the 50 Hz interference component.

The filter parameters were:

| Parameter | Value |
|---|---|
| Center frequency | 50 Hz |
| Bandwidth | 5 Hz |
| G0 | 1 |
| Gc | 0.1 |

![Figure 2: ECG time-domain and frequency-domain comparison](ecg_notch.png)

*Figure 2. Comparison of the clean, noisy and filtered ECG-like signals in the time domain, together with the frequency spectra before and after filtering.*

The time-domain plot shows that the 50 Hz oscillations are substantially reduced after filtering. The filtered waveform more closely follows the clean ECG-like reference, although some baseline wander and residual noise remain.

The frequency-domain plot confirms that the component at 50 Hz is attenuated. This demonstrates the notch filter's ability to suppress a selected interference frequency while retaining much of the underlying waveform.

The measured mean squared error (MSE) was reduced by approximately 87.7% in this synthetic experiment.

### Tone-like Signal — Resonance Filtering

A resonance filter was applied to amplify the desired 100 Hz component of a synthetic audio-like signal.

The filter parameters were:

| Parameter | Value |
|---|---|
| Center frequency | 100 Hz |
| Bandwidth | 10 Hz |
| G0 | 1 |
| Gc | 5 |

![Figure 3: Tone time-domain and frequency-domain comparison](tone_resonance.png)

*Figure 3. Comparison of the clean, noisy and filtered tone-like signals in the time domain, together with the frequency spectra before and after resonance filtering.*

The frequency-domain plot shows a significant increase in the amplitude around 100 Hz after filtering. The interference components at 40 Hz, 180 Hz and 300 Hz remain visible because the resonance filter is designed to amplify a selected frequency region rather than remove all unwanted components.

In the time domain, the filtered signal does not necessarily resemble the clean reference. This is expected because resonance filtering amplifies the desired frequency component without completely suppressing the other components.

The numerical validation showed that the actual gain at 100 Hz was approximately 8.014 instead of the requested gain of 5, corresponding to an error of +60.29%.

This is an important limitation of the approximate pole-zero placement method for the selected parameters.

### Influence of Bandwidth

The frequency responses were also compared using two bandwidth values for each filter type.

![Figure 4: Comparison of notch and resonance filter bandwidths](bandwidth_comparison.png)

*Figure 4. Comparison of frequency responses for notch filters with bandwidths of 5 Hz and 0.5 Hz, and resonance filters with bandwidths of 10 Hz and 1 Hz.*

For the notch filters, reducing the bandwidth from 5 Hz to 0.5 Hz produces a narrower attenuation region around 50 Hz. The response outside this region remains closer to the desired baseline gain.

For the resonance filters, reducing the bandwidth from 10 Hz to 1 Hz produces a narrower peak around 100 Hz. The actual center gain also becomes closer to the requested value, decreasing from approximately 8.014 to 5.247.

However, the relative bandwidth error increases from approximately 1.52% to 3.90% in the resonance examples. Therefore, narrower bandwidth does not necessarily improve the accuracy of every filter characteristic.

### Discussion

The experiments demonstrate that the implemented filters perform their intended qualitative functions: the notch filter attenuates a selected frequency, while the resonance filter amplifies a selected frequency region.

All four tested configurations were stable, but their gain accuracy varied considerably. The results highlight the importance of evaluating the actual frequency response instead of assuming that the approximate design equations reproduce all requested specifications exactly.

The experiments also illustrate the difference between selective attenuation and selective amplification. A resonance filter should not be interpreted as a general noise-removal filter.

## 6. Limitations

- The design equations are approximations, and the errors depend on the requested parameters.
- Only four filter configurations were evaluated numerically.
- The signal examples are **synthetic, not real-world recordings**. If the assignment requires real-life signals, recordings or externally sourced measured data must be added before submission.
- The validation measures frequency response numerically; its resolution and interpolation affect the reported precision.
- The notebook's MSE uses a known synthetic reference, which would not ordinarily be available for a real measurement.

## 7. Reproducing the experiments

Requires Python with `numpy`, `scipy`, `matplotlib`, and Jupyter (`notebook` or `jupyterlab`). Run from the project directory:

```bash
python generate_signals.py
python validate_filter.py
jupyter notebook filter_analysis.ipynb
```

In the notebook, select **Restart Kernel and Run All Cells**. Generated WAV files are stored in `signals/`; filtered WAV outputs are written to `output/`. Generated data and outputs can be excluded from version control because they are reproducible.

## Reference

Van Waterschoot, T., and Moonen, M. (2007). *A Pole-Zero Placement Technique for Designing Second-Order IIR Parametric Equalizer Filters.* Refer to the paper supplied with the assignment for the complete derivation and publication details.
