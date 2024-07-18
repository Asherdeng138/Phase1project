clear all;
close all;
clc;

N = 10^4; % number of bits or symbols
Eb_N0_dB = -3:35; % multiple Eb/N0 values

% BPSK Simulation
nErr_BPSK = zeros(size(Eb_N0_dB));
ip = rand(1,N)>0.5; % generating 0,1 with equal probability
s_BPSK = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 1

% Generate s_QPSK for the entire simulation beforehand
s_QPSK = qammod(double(ip),4,'UnitAveragePower',true); % QPSK modulation, produces N/2 symbols

for ii = 1:length(Eb_N0_dB)
    n_BPSK = 1/sqrt(2)*(randn(1,N) + 1j*randn(1,N)); % white gaussian noise, 0dB variance 
    h_BPSK = 1/sqrt(2)*(randn(1,N) + 1j*randn(1,N)); % Rayleigh channel
    y_BPSK = h_BPSK.*s_BPSK + 10^(-Eb_N0_dB(ii)/20)*n_BPSK; % Channel and noise addition
    yHat_BPSK = y_BPSK./h_BPSK; % equalization
    ipHat_BPSK = real(yHat_BPSK)>0; % receiver - hard decision decoding
    nErr_BPSK(ii) = sum(ip ~= ipHat_BPSK); % counting the errors
end
simBer_BPSK = nErr_BPSK/N; % simulated BER

% QPSK Simulation
nErr_QPSK = zeros(size(Eb_N0_dB));

% 修正: 确保s_QPSK生成时与h_QPSK和n_QPSK匹配
% 注意: 假设N是偶数
s_QPSK = qammod(double(ip(1:2:end)*2 + ip(2:2:end)), 4, 'UnitAveragePower', true); % QPSK modulation

for ii = 1:length(Eb_N0_dB)
    % Adjust noise and channel size for QPSK
    n_QPSK = 1/sqrt(2)*(randn(1,N/2) + 1j*randn(1,N/2)); % AWGN for QPSK
    h_QPSK = 1/sqrt(2)*(randn(1,N/2) + 1j*randn(1,N/2)); % Rayleigh channel for QPSK
    
    % Correctly sized operation for channel and noise addition for QPSK
    y_QPSK = h_QPSK .* s_QPSK + 10^(-Eb_N0_dB(ii)/20) * n_QPSK;
    
    % Equalization for QPSK
    yHat_QPSK = y_QPSK ./ h_QPSK;
    
    % QPSK Demodulation
    ipHat_QPSK = qamdemod(yHat_QPSK, 4, 'UnitAveragePower', true);

    % Convert demodulated symbols back to bits
    decodedBits_QPSK = de2bi(ipHat_QPSK, 'left-msb', 2); 
    decodedBits_QPSK = reshape(decodedBits_QPSK', 1, []); 
    
    % Counting the errors for QPSK, adjusted for the total number of bits
    nErr_QPSK(ii) = sum(ip(1:N) ~= decodedBits_QPSK(1:N));
end
simBer_QPSK = nErr_QPSK / N; % average BER for QPSK


% Theoretical BER
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer_BPSK = 0.5*(1-sqrt(EbN0Lin./(EbN0Lin+1))); % Theoretical BER for BPSK over Rayleigh
theoryBer_QPSK = 0.5*(1-sqrt((EbN0Lin/2)./(EbN0Lin/2+1))); % Corrected Theoretical BER for QPSK over Rayleigh


% Plot
figure;
semilogy(Eb_N0_dB, theoryBer_BPSK, 'bp-', 'LineWidth', 2, 'DisplayName', 'BPSK Rayleigh Theory');
hold on;
semilogy(Eb_N0_dB, simBer_BPSK, 'mx-', 'LineWidth', 2, 'DisplayName', 'BPSK Rayleigh Simulation');
semilogy(Eb_N0_dB, theoryBer_QPSK, 'gp-', 'LineWidth', 2, 'DisplayName', 'QPSK Rayleigh Theory');
semilogy(Eb_N0_dB, simBer_QPSK, 'rx-', 'LineWidth', 2, 'DisplayName', 'QPSK Rayleigh Simulation');

legend('show');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
title('BER for BPSK and QPSK Modulation in Rayleigh Channel');
axis([-3 35 10^-5 0.5]);
grid on;
