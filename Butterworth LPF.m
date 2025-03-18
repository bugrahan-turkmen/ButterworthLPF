% DSP

% Sample script
% Step 1: Sampling

% Parameters 
Fs = 16000; % Sampling frequency 
T = 1 / Fs; % Sampling period
% Amplitudes
A1 = 1; 
A2 =0.7; 

% Frequencies
F1 =300;
F2 =1100; 

% Phases
phi1 =1;
phi2 = 1.5; 

% Time vector for 2048 samples
n = 0:2047;
t = n * T; % Continuous time equivalent

% Generate DT signal x[n]
xn = A1*cos(2*(pi)*F1*t+phi1) + A2*cos(2*(pi)*F2*t+phi2);

% Plot x[n] for n between 256 and 319. Length=64
figure;
n1=256:319;
stem(n1, xn(n1));
xlabel('n');
ylabel('Amplitude');
title('Discrete-Time Signal x[n]');
grid on;

% Step 2: Fourier Transform of x[n]

% Compute FFT
Xk = fft(xn,1024);
N = length(Xk);
frequencies = (0:N-1) * (Fs / N);

% Magnitude and phase
magnitude = abs(Xk);
phase = angle(Xk);

% Plot magnitude spectrum
figure;
plot(frequencies, magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum');
grid on;

% Plot phase spectrum
figure;
plot(frequencies, phase);
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Spectrum');
grid on;

% Determine peaks.
[~, indices] = sort(magnitude(1:N/2) , 'descend');
F1_peak = frequencies(indices(2));
F2_peak = frequencies(indices(3));

% Percentage errors
error_F1 = abs((F1 - F1_peak) / F1) * 100;
error_F2 = abs((F2 - F2_peak) / F2) * 100;

% Step 3: Low-Pass Filter Design

% Design a low-pass filter (adjust order and cutoff frequency as needed)
filter_order = 8;

% Cutoff frequency= Fc/(Fs/2). A low-pass filter with cutoff (in Hz) 𝐹𝑐 = 0.1 × 𝐹𝑠/2.

%Fc = 400;
Fc = 0.05 * (Fs/2);
cutoff_frequency = Fc/(Fs/2);

[b, a] = butter(filter_order, cutoff_frequency, 'low');


% Frequency response
[H, W] = freqz(b, a, 1024, Fs);

% Plot frequency response
figure;
subplot(2, 1, 1);
plot(W, 20 * log10(abs(H)));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Response');
grid on;
subplot(2, 1, 2);
plot(W, angle(H));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Response');
grid on;

% Find -3dB cutoff frequency

cutoff_idx = find(20*log10(abs(H)) <= -3, 1);
cutoff_freq = W(cutoff_idx);
fprintf('-3dB Cutoff Frequency: %.2f Hz\n', cutoff_freq);

% Step 4: System Output

% Filter the input signal
yn = filter(b, a, xn);



% Plot output signal
figure;
stem(n1, yn(n1));
xlabel('n');
ylabel('Amplitude');
title(' Output (y[n]) Signals');
grid on;

% FFT of the output
Yk = fft(yn, 1024);
magnitude_y = abs(Yk);
phase_y = angle(Yk);

% Plot magnitude spectrum of output
figure;
plot(frequencies, 20 * log10(magnitude_y));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Spectrum of Output y[n]');
grid on;

% Plot phase spectrum of output
figure;
plot(frequencies, phase_y);
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Spectrum of Output y[n]');
grid on;

% Determine index for the second sinusoid
ik2 = round(N*F2/Fs);

% Determine frequency specified by the index
output_freq_F2 = frequencies(ik2);

% Determine the magnitude in dB at the given frequency
output_magnitude_F2 = 20 * log10(magnitude_y(ik2));
fprintf('Second Sinusoid Frequency: %.2f Hz\n', output_freq_F2);
fprintf('Second Sinusoid Magnitude: %.2f dB\n', output_magnitude_F2);

% Adjust filter if necessary
if output_magnitude_F2 > -6
 disp('Adjusting filter order or cutoff frequency to meet specification.');
end