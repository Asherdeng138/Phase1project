function Pe = symbol_error_prob(M, EbN0_dB)
    % Convert SNR from dB to linear scale
    gamma_b = 10.^(EbN0_dB/10);
    
    % MATLAB's Q-function equivalent
    Q = @(x) 0.5 * erfc(x / sqrt(2));
    
    % BPSK
    if M == 2
        Pe = Q(sqrt(2 * gamma_b));
    % QPSK
    elseif M == 4
        Pe = 2 * Q(sqrt(2 * gamma_b)) - (Q(sqrt(2 * gamma_b)))^2;
    % M-ary PSK for M > 2
    else
        Pe = 2 * Q(sqrt(2 * log2(M) * gamma_b * sin(pi/M)));
    end
end


% Define the SNR range and modulation orders
snr_dB = linspace(-4, 24, 1000);
M_values = [2, 4];

% Plot the symbol error probability
figure;
for M = M_values
    ser = arrayfun(@(x) symbol_error_prob(M, x), snr_dB);
    semilogy(snr_dB, ser, 'DisplayName', sprintf('M = %d', M));
    hold on;
end

title('Probability of Symbol Error for PSK Signals');
xlabel('SNR per bit, \gamma_b (dB)');
ylabel('Probability of Symbol Error, P_e');
legend('show');
grid on;
ylim([1e-5, 1]);
