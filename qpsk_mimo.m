clear all;
close all;
clc;

% Initialize parameters
SymbolMapping_ = 'Gray'; % Symbol mapping method
SNR_dB = -5:1:15; % SNR range in dB
bit_number = 10^6; % Number of bits
Frame = bit_number / 4; % Number of frames

% Create QPSK modulator and demodulator objects
qpsk_mo = comm.QPSKModulator('SymbolMapping', SymbolMapping_, 'BitInput', true);
qpsk_demo = comm.QPSKDemodulator('SymbolMapping', SymbolMapping_, 'BitOutput', true);

% Initialize BER storage
BER_MISO = zeros(1, length(SNR_dB));
BER_SIMO = zeros(1, length(SNR_dB));
BER_MIMO = zeros(1, length(SNR_dB));

tic;

for i = 1:length(SNR_dB)
    ErrorNum_MISO = 0;
    ErrorNum_SIMO = 0;
    ErrorNum_MIMO = 0;

    for j = 1:Frame
        % Generate random bits
        bit_in = randi([0 1], 1, 4).';

        % Modulate bits
        bit_channel = qpsk_mo(bit_in);
        x1 = bit_channel(1);
        x2 = bit_channel(2);

        % --------------------------- MISO Simulation ---------------------------
        % Space-time coding
        X_1 = [x1 -conj(x2)];
        X_2 = [x2 conj(x1)];
        % Pass through AWGN channel
        X_1_channel = awgn(X_1, SNR_dB(i));
        X_2_channel = awgn(X_2, SNR_dB(i));
        % Receive signal
        R_MISO = X_1_channel + X_2_channel;
        % Space-time decoding
        X1_MISO = (R_MISO(1) + conj(R_MISO(2))) / 2;
        X2_MISO = (R_MISO(1) - conj(R_MISO(2))) / 2;
        bit_channel_awgn_MISO = [X1_MISO X2_MISO].';

        % Demodulate
        bit_out_MISO = qpsk_demo(bit_channel_awgn_MISO);
        % Calculate the number of bit errors
        ErrorNum_MISO = ErrorNum_MISO + sum(bit_in ~= bit_out_MISO);
        
        % --------------------------- SIMO Simulation ---------------------------
        % Generate 2 different Rayleigh channels
        h = (randn(2, 1) + 1j*randn(2, 1)) / sqrt(2);
        % Pass through AWGN channel
        R_SIMO = zeros(2, 2);
        for k = 1:2
            R_SIMO(k, :) = awgn(h(k) * [x1; x2].', SNR_dB(i), 'measured');
        end
        % Combine received signals
        R_combined = sum(conj(h) .* R_SIMO, 1) / sum(abs(h).^2);
        % Space-time decoding
        X_SIMO = R_combined.';
        bit_channel_awgn_SIMO = X_SIMO;

        % Demodulate
        bit_out_SIMO = qpsk_demo(bit_channel_awgn_SIMO);
        % Calculate the number of bit errors
        ErrorNum_SIMO = ErrorNum_SIMO + sum(bit_in ~= bit_out_SIMO);

        % --------------------------- MIMO Simulation ---------------------------
        % Pass through AWGN channel
        X1_1_channel = awgn(X_1, SNR_dB(i));
        X1_2_channel = awgn(X_2, SNR_dB(i));
        % Receive signals
        R1_MIMO = X_1_channel + X_2_channel;
        R2_MIMO = X1_1_channel + X1_2_channel;
        % Space-time decoding
        X11_MIMO = (R1_MIMO(1) + conj(R1_MIMO(2))) / 2;
        X21_MIMO = (R1_MIMO(1) - conj(R1_MIMO(2))) / 2;
        bit_channel_awgn1_MIMO = [X11_MIMO X21_MIMO].';

        X12_MIMO = (R2_MIMO(1) + conj(R2_MIMO(2))) / 2;
        X22_MIMO = (R2_MIMO(1) - conj(R2_MIMO(2))) / 2;
        bit_channel_awgn2_MIMO = [X12_MIMO X22_MIMO].';

        bit_channel_awgn_MIMO = (bit_channel_awgn1_MIMO + bit_channel_awgn2_MIMO) / 2;

        % Demodulate
        bit_out_MIMO = qpsk_demo(bit_channel_awgn_MIMO);
        % Calculate the number of bit errors
        ErrorNum_MIMO = ErrorNum_MIMO + sum(bit_in ~= bit_out_MIMO);
    end

    % Calculate BER
    BER_MISO(i) = ErrorNum_MISO / (Frame * 4);
    BER_SIMO(i) = ErrorNum_SIMO / (Frame * 4);
    BER_MIMO(i) = ErrorNum_MIMO / (Frame * 4);

    % Display progress
    fprintf('SNR = %d dB, MISO BER = %.6f, SIMO BER = %.6f, MIMO BER = %.6f\n', SNR_dB(i), BER_MISO(i), BER_SIMO(i), BER_MIMO(i));
end

toc;

% Plot results
figure;
semilogy(SNR_dB, BER_MISO, 'b-o', 'DisplayName', 'MISO-QPSK Practical');
hold on;
semilogy(SNR_dB, BER_SIMO, 'g-*', 'DisplayName', 'SIMO-QPSK Practical');
semilogy(SNR_dB, BER_MIMO, 'r-s', 'DisplayName', 'MIMO-QPSK Practical');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance of MISO, SIMO, and MIMO QPSK over AWGN Channel');
legend('show');
grid on;
