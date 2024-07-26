clear all; close all; clc;

% Initialize parameters
N = 10^6; % Number of bits
M_BPSK = 2;
M_QPSK = 4; % Modulation order for QPSK
snr_db = 0:1:10; % SNR range in dB

% Generate BPSK and QPSK data
binary_data = randsrc(1, N, [0, 1]);
BPSK_data = 2*binary_data - 1;
QPSK_data = qammod(binary_data(1:2:end)*2 + binary_data(2:2:end), M_QPSK, 'UnitAveragePower', true);

% Preallocate BER arrays for practical and theoretical values
BER_BPSK_prac = zeros(1, length(snr_db));
BER_QPSK_prac = zeros(1, length(snr_db));
BER_BPSK_theo = qfunc(sqrt(2*10.^(snr_db/10)));
BER_QPSK_theo = qfunc(sqrt(10.^(snr_db/10)));

% Calculate BER for each SNR value
for k = 1:length(snr_db)
    SNR_linear = 10.^(snr_db(k)/10);
    noise_power = 1 ./ SNR_linear;
    noise_std = sqrt(noise_power / 2);
    noise = noise_std .* (randn(1, N) + 1i*randn(1, N));
    
    % BPSK
    y_BPSK = BPSK_data + noise;
    decoded_BPSK = real(y_BPSK) > 0;
    BER_BPSK_prac(k) = mean(decoded_BPSK ~= binary_data);
    
    % QPSK
    y_QPSK = QPSK_data + noise(1:N/2); % Apply noise to each QPSK symbol
    decoded_QPSK = qamdemod(y_QPSK, M_QPSK, 'UnitAveragePower', true);
    decoded_bits = reshape(de2bi(decoded_QPSK, log2(M_QPSK), 'left-msb')', 1, []);
    BER_QPSK_prac(k) = mean(decoded_bits(1:N) ~= binary_data(1:N)); % Compare the first N bits
end

% Plot the results
figure;
semilogy(snr_db, BER_BPSK_prac, 'b-o', 'DisplayName', 'BPSK Practical');
hold on;
semilogy(snr_db, BER_QPSK_prac, 'r-s', 'DisplayName', 'QPSK Practical');
semilogy(snr_db, BER_BPSK_theo, 'b--', 'DisplayName', 'BPSK Theoretical');
semilogy(snr_db, BER_QPSK_theo, 'r--', 'DisplayName', 'QPSK Theoretical');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance of BPSK and QPSK over AWGN Channel');
legend('show');
grid on;

