A = load('probs_alpha.txt');
B = load('probs_richardson.txt');
x = 1: 1: size(B, 1);
plot(x, A(:, 2), 'b.');
grid on;
hold on;
plot(x, B(:, 2), 'r.');
title('Probability distribution with different processes over 70 cell Yung diagrams');
legend('Alpha = 0.3', 'Richardson');
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'probs.jpg');

A = load('probs_alpha.txt');
x = 1: 1: size(A, 1);
plot(x, sort(A(:, 2)), 'b.');
hold on;
grid on;
C = load('probs_beta.txt');
x = 1: 1: size(C, 1);
plot(x, sort(C(:, 2)), 'r.');
legend('Alpha = 0.16', 'Вероятность обратно пропорциональна числу входящих дуг');
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'InverseProbs.jpg');


figure;
C = load('Freqs_400_1000buckets_alpha.txt');
%C = sort(C);
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


% assymptotic shape
A = load('100000Cells_Richardson.txt');
x = 1: 1: size(A, 2);
bar(x, A, 'r');
grid on;
hold on;
A = load('100000Cells_ALPHA.txt');
x = 1: 1: size(A, 2);
bar(x, A, 'b');
grid on;
hold on;
A = load('100000Cells_PLANSHEREL_GAMMA.txt');
x = 1: 1: size(A, 2);
bar(x, A, 'b');
grid on;
hold on;


A = load('1000000Cells_BETA.txt');
x = 1: 1: size(A, 2);
plot(x, A, 'g');