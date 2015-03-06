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
A = load('1000000Cells_Richardson.txt');
x = 1: 1: size(A, 2);
plot(x, A, 'b');
grid on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'Richardson_assympt.jpg');
hold on;

A = load('1000000Cells_ALPHA.txt');
x = 1: 1: size(A, 2);
bar(x, A, 'b');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients estimation
A = load('EstimationCoefs.txt');
c = A(1, :);
X = A(2:end,:);
y = X * c';
x = 1:1:size(X,1);
plot(x, y, '*b');
min(y)
grid on;
hold on;
B = load('EstimationDistances.txt');
Y = B(:, 2);
plot(x, Y, '*r');
legend('Estimation', 'Real distances');
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'MetricEstimation20.jpg');
hold off;

x = 1:1:size(c, 2);
plot(x, c, '*b');
grid on;
legend('Оценка c(i)');
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'CoefficientsEstimation20.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 3d assymptotic shape
A = importdata('3dHooksRandomFast.txt');
A(isnan(A)) = 0;
figure;
%A = A .* 100;
h = bar3(A);
[r c] = size(A);
for i = 1:c
    zdata = [];
    for j = 1:r
        zdata = [zdata; ones(6,4)*A(j,i)];
    end
    set(h(i),'Cdata',zdata)
end
colormap jet
for i = 1:numel(h)
  index = logical(kron(A(:,i) == 0,ones(6,1)));
  zData = get(h(i),'ZData');
  zData(index,:) = nan;
  set(h(i),'ZData',zData);
end
colorbar
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, '3DHooks100000.jpg');


x = ones(size(A, 1) * size(A, 2), 1);
y = ones(size(A, 1) * size(A, 2), 1);
z = ones(size(A, 1) * size(A, 2), 1);
ind = 1;
for i = 1 : 1 : size(A, 1)
    for j = 1: 1 : size(A, 2);
        if (A(i, j) ~= 0)
            x(ind) = i;
            y(ind) = j;
            z(ind) = A(i, j);
            ind = ind + 1;
        end
    end
end
ind = ind - 1;
x = x(1:ind);
y = y(1:ind);
z = z(1:ind);

figure;
F = TriScatteredInterp(x, y, z);
ti = 0:0.1:135;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);
mesh(qx,qy,qz);

t = load('3dtime.txt');
x = t(:, 1);
y = t(:, 2);
coeff = polyfit(x, y, 2);
xfit = linspace(min(x),max(x),20);
yfit = polyval(coeff,xfit);
h = plot(x,y,'*',xfit,yfit,'r');
set(h,'LineWidth',3);
xlabel('Number of cells');
ylabel('Time, seconds');
grid on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'AssymptHooksTime.jpg');

polyval(coeff, 100000);