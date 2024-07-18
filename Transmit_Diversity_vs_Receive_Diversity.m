% MATLAB Simulation for Transmit vs. Receive Diversity
clear all;
close all;
clc;

% Simulation Parameters
frmLen = 100;       % Frame length (number of symbols per frame)
numPackets = 1000;  % Number of packets
EbNo = 0:2:20;      % Eb/No values in dB
N = 2;              % Number of Tx antennas
M = 2;              % Number of Rx antennas
P = 2;              % Modulation order (BPSK)

% Create OSTBC Encoder and Combiner System objects
ostbcEnc = comm.OSTBCEncoder;
ostbcComb = comm.OSTBCCombiner;

% Convert Eb/No values to SNR values (unit power signals for BPSK)
SNR = convertSNR(EbNo, "ebno", "BitsPerSymbol", 1);

% Create ErrorRate calculator System objects to evaluate BER
errorCalc1 = comm.ErrorRate;
errorCalc2 = comm.ErrorRate;
errorCalc3 = comm.ErrorRate;

% Set up random stream for repeatability
s = rng(55408);

% Pre-allocate variables for speed
H = zeros(frmLen, N, M);
ber_noDiver = zeros(3, length(EbNo));
ber_Alamouti = zeros(3, length(EbNo));
ber_MaxRatio = zeros(3, length(EbNo));
ber_thy2 = zeros(1, length(EbNo));

% Set up a figure for visualizing BER results
fig = figure;
grid on;
ax = fig.CurrentAxes;
hold(ax, 'on');
ax.YScale = 'log';
xlim(ax, [EbNo(1), EbNo(end)]);
ylim(ax, [1e-4 1]);
xlabel(ax, 'Eb/No (dB)');
ylabel(ax, 'BER');
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
fig.Name = 'Transmit vs. Receive Diversity';
title(ax, 'Transmit vs. Receive Diversity');
set(fig, 'DefaultLegendAutoUpdate', 'off');
fig.Position = figposition([15 50 25 30]);

% Loop over several Eb/No points
for idx = 1:length(EbNo)
    reset(errorCalc1);
    reset(errorCalc2);
    reset(errorCalc3);
    
    % Loop over the number of packets
    for packetIdx = 1:numPackets
        % Generate data vector per frame
        data = randi([0 P-1], frmLen, 1);

        % Modulate data using BPSK
        modData = pskmod(data, P);

        % Alamouti Space-Time Block Encoder
        encData = ostbcEnc(modData);

        % Create the Rayleigh distributed channel response matrix for 2x2 MIMO
        H(1:N:end, :, :) = (randn(frmLen/2, N, M) + 1i*randn(frmLen/2, N, M)) / sqrt(2);
        H(2:N:end, :, :) = H(1:N:end, :, :);  % Assume constant for 2 symbol periods

        % Extract part of H to represent the 1x1, 2x1 and 1x2 channels
        H11 = H(:, 1, 1);
        H21 = H(:, :, 1) / sqrt(2);
        H12 = squeeze(H(:, 1, :));

        % Pass through the channels
        chanOut11 = H11 .* modData;
        chanOut21 = sum(H21 .* encData, 2);
        chanOut12 = H12 .* repmat(modData, 1, 2);

        % Add AWGN to the channel outputs
        rxSig11 = awgn(chanOut11, SNR(idx));
        rxSig21 = awgn(chanOut21, SNR(idx));
        rxSig12 = awgn(chanOut12, SNR(idx));

        % Alamouti Space-Time Block Combiner
        decData = ostbcComb(rxSig21, H21);

        % ML Detector (minimum Euclidean distance)
        demod11 = pskdemod(rxSig11 .* conj(H11), P);
        demod21 = pskdemod(decData, P);
        demod12 = pskdemod(sum(rxSig12 .* conj(H12), 2), P);

        % Calculate and update BER for current Eb/No value
        ber_noDiver(:, idx) = errorCalc1(data, demod11);
        ber_Alamouti(:, idx) = errorCalc2(data, demod21);
        ber_MaxRatio(:, idx) = errorCalc3(data, demod12);
    end

    % Calculate theoretical second-order diversity BER for current Eb/No
    ber_thy2(idx) = berfading(EbNo(idx), 'psk', 2, 2);

    % Plot results
    semilogy(ax, EbNo(1:idx), ber_noDiver(1, 1:idx), 'r*', ...
             EbNo(1:idx), ber_Alamouti(1, 1:idx), 'go', ...
             EbNo(1:idx), ber_MaxRatio(1, 1:idx), 'bs', ...
             EbNo(1:idx), ber_thy2(1:idx), 'm');
    legend(ax, 'No Diversity (1Tx, 1Rx)', 'Alamouti (2Tx, 1Rx)',...
           'Maximal-Ratio Combining (1Tx, 2Rx)', 'Theoretical 2nd-Order Diversity');

    drawnow;
end

% Perform curve fitting and replot the results
fitBER11 = berfit(EbNo, ber_noDiver(1, :));
fitBER21 = berfit(EbNo, ber_Alamouti(1, :));
fitBER12 = berfit(EbNo, ber_MaxRatio(1, :));
semilogy(ax, EbNo, fitBER11, 'r', EbNo, fitBER21, 'g', EbNo, fitBER12, 'b');
hold(ax, 'off');

% Restore default random stream
rng(s);
