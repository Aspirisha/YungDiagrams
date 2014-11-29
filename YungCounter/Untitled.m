A = load('probs_alpha.txt');
B = load('probs_richardson.txt');
x = 1: 1: size(A, 1);
plot(x, A(:, 2), 'b.');
grid on;
hold on;
plot(x, B(:, 2), 'r.');
title('Probability distribution with different processes over 70 cell Yung diagrams');
legend('Alpha = 0.3', 'Richardson');
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'probs.jpg');

figure;
C = load('Freqs_400_1000buckets_alpha.txt');
C = sort(C);
x = 1: 1: size(C);
bar(x, C, 'b');
title('Frequencies for 1000 diagram buckets for diagrams with 400 cells.');
legend('Alpha = 0.16');
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'Freqs_400_1000buckets_alpha.jpg');
grid on;

figure;
C = load('Freqs_400_1000buckets_richardson.txt');
C = sort(C);
x = 1: 1: size(C);
bar(x, C, 'b');
title('Frequencies for 1000 diagram buckets for diagrams with 400 cells.');
legend('Richardson');
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'Freqs_400_1000buckets_richardson.jpg');
grid on;