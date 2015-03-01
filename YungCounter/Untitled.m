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
legend('Alpha = 0.16', '����������� ������� ��������������� ����� �������� ���');
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
A = load('1000000Cells_Richardson.txt');
x = 1: 1: size(A, 2);
plot(x, A, 'b');
grid on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'Richardson_assympt.jpg');
hold on;

A = load('1000000Cells_ALPHA.txt');
x = 1: 1: size(A, 2);
plot(x, A, 'r');
grid on;
hold on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'ALPHA_assympt.jpg');

A = load('1000000Cells_PLANSHEREL.txt');
x = 1: 1: size(A, 2);
plot(x, A, 'g');
grid on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'Plansherel_assympt.jpg');
hold on;


A = load('1000000Cells_BETA.txt');
x = 1: 1: size(A, 2);
bar(x, A, 'b');
grid on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'Beta_assympt.jpg');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simplex distances
D = load('../KantorovichDistanceCounter/SimplexDistances.txt');
x = D(:, 1);
y = D(:, 2);
z = D(:, 3);

figure;
F = TriScatteredInterp(x, y, z);
ti = -2:0.01:2;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);
mesh(qx,qy,qz);
%set(gcf, 'InvertHardCopy', 'off');
%saveas(gcf, 'SimplexDists1.jpg');
hold on;

x = x - 2;
plot3(x,y,z,'o');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kantorovich ball
A = load('KantorovichBall.txt');
plot(A(:,1), A(:,2), '-r');

A = load('AllDistances20.txt');
plot(A(:,1), A(:,2), '-r');
grid on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'DiagramsDists20.jpg');
