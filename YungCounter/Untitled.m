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