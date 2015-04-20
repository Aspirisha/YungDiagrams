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
A = load('BallRandomPlansherel20_Radius_0_1.txt');
plot(A(:,1), A(:,2), '-r');

A = load('AllDistances20FromLine.txt');
plot(A(:,1), A(:,2), '-r');
grid on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'DiagramsDists20FromLine.jpg');

A = load('AllDistances20.txt');
plot(A(:,1), A(:,2), '-r');
grid on;
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'DiagramsDists20.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients estimation
figure;
A = load('EstimationCoefsFabs_35.txt');
c = A(1, :);
X = A(2:end,:);
y = X * c';
x = 1:1:size(X,1);
plot(x, y, '*b');
min(y)
grid on;
hold on;
B = load('EstimationDistancesFabs_35.txt');
Y = B(:, 2);
plot(x, Y, '*r');
legend('Estimation', 'Real distances');
%set(gcf, 'InvertHardCopy', 'off');
%saveas(gcf, 'MetricEstimation20Fabs.jpg');
hold off;

figure;
x = 1:1:size(c, 2);
plot(x, c, '*b');
grid on;
legend('Оценка c(i)');
%set(gcf, 'InvertHardCopy', 'off');
%saveas(gcf, 'CoefficientsEstimation20Fabs.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 3d assymptotic shape
A = importdata('3dHooksRandomFast1000000.txt');
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

% scaled in-corner
A = importdata('3dHooksRandomFast1000000.txt');
A(isnan(A)) = 0;
ind = 1;
s = 1000000 ^ (-1/3);
clearvars x y z;
for i = 1 : 1: size(A, 1)
    for j = 1 : 1: size(A, 2)
        if (A(i,j) ~= 0)
            x(ind) = i * s;
            y(ind) = j * s;
            z(ind) = A(i, j) * s;
            ind = ind + 1;
        end
    end;
end
plot3(x,y,z,'r*');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D hooks asymptotic function

% scaled and rotated
A = importdata('3dHooksRandomFast1000000.txt');
A(isnan(A)) = 0;
ind = 1;
s = 1000000 ^ (-1/3);
clearvars x y z;
for i = 1 : 1: size(A, 1)
    for j = 1 : 1: size(A, 2)
        x1 = i * s;
        y1 = j * s;
        z1 = A(i, j) * s;
        
        if (z1 ~= 0)
            x(ind) = (x1 - y1) / sqrt(2);
            y(ind) = (2 * z1 - x1 - y1) / sqrt(6);
            z(ind) = (x1 + y1 + z1) / sqrt(3);
            ind = ind + 1;
        end
    end;
end
plot3(x,y,z,'r*');
grid on;

min(z) - 9/(pi*pi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time for 3D hooks
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

polyval(coeff, 1000000)




% integration
syms t

int(((2/3)*asin((1/2)*t)*(t*asin((1/2)*t)+sqrt(-t^2+4))/pi^2-t)/(t^2+(1/3)*(2*(t*asin((1/2)*t)+sqrt(-t^2+4))/pi)^2), t)

fun = @(t) ((2/3)*asin(0.5.*t).*(t.*asin(0.5*t)+sqrt(-t.^2+4))/pi^2-t)./(t.^2+(1/3)*(2.*(t.*asin(0.5*t)+sqrt(-t.^2+4))/pi).^2)
q = quadgk(@(t)fun(t).*fun(t),-2,2)
%
[q, n] = quad('((2/3)*asin(0.5.*t).*(t.*asin(0.5*t)+sqrt(-t.^2+4))/pi^2-t)./(t.^2+(1/3)*(2.*(t.*asin(0.5*t)+sqrt(-t.^2+4))/pi).^2)', -2, 2)


x0 = 0;
y0 = 0;
f = @(t) 2/pi * (t.*asin(t/2) + sqrt(4 - t.*t));
df = @(t) 2/pi * asin(t/2);
x = @(t, fi) t.*cos(fi) - f(t).*sin(fi)/sqrt(3);
y = @(t, fi) -(t.*sin(fi) + f(t).*cos(fi)/sqrt(3));
z = @(t) f(t) * sqrt(2/3);
dx = @(t, fi) cos(fi) - df(t).*sin(fi)/sqrt(3);
dy = @(t, fi) -(sin(fi) + df(t).*cos(fi)/sqrt(3));
F = @(t, x0, y0, fi) 1/(2*pi)*z(t).*(-dx(t, fi).*(y(t, fi) - y0) + dy(t, fi).*(x(t, fi) - x0))./((x(t, fi)-x0).^2 + (y(t, fi)-y0).^2);

clearvars X Y Z
ind = 1
for i = -2:0.01:2
    X(ind) = x(i, pi/4);
    Y(ind) = y(i, pi/4);
    Z(ind) = z(i);
    ind = ind + 1
end

for i = -2:0.01:2
    X(ind) = x(i, -5*pi/12);
    Y(ind) = y(i, -5*pi/12);
    Z(ind) = z(i);
    ind = ind + 1
end

for i = -2:0.01:2
    X(ind) = x(i, -13*pi/12);
    Y(ind) = y(i, -13*pi/12);
    Z(ind) = z(i);
    ind = ind + 1
end

plot3(X,Y,Z)
grid on;
hold on;

% scaled and rotated
A = importdata('3dHooksRandomFast1000000.txt');
A(isnan(A)) = 0;
ind = 1;
s = 1000000 ^ (-1/3);
clearvars x y z;
clearvars X Y Z;
for i = 1 : 1: size(A, 1)
    for j = 1 : 1: size(A, 2)
        x1 = i * s;
        y1 = j * s;
        z1 = A(i, j) * s;
        
        if (z1 ~= 0)
            x(ind) = (x1 - y1) / sqrt(2);
            y(ind) = (2 * z1 - x1 - y1) / sqrt(6);
            z(ind) = (x1 + y1 + z1) / sqrt(3);
            X(ind) = (x(ind) + y(ind)) / sqrt(2);
            Y(ind) = (-x(ind) + y(ind)) / sqrt(2);
            Z(ind) = z(ind);
            ind = ind + 1;
        end
        
    end;
end
plot3(X,Y,Z,'r*');
grid on;
hold on

g = @(t) F(t,0,0,pi/4); 
q = quad(@(t)g(t),-2,2);
g = @(t) F(t,0,0,-5*pi/12); 
q = q + quad(@(t)g(t),-2,2);
g = @(t) F(t,0,0,-13*pi/12); 
q = q + quad(@(t)g(t),-2,2);
q

n = size(X, 2)
ind = 1;
Zint = zeros(1, n);
Xint = zeros(1, n);
Yint = zeros(1, n);
for j = 1:1:n
    X0 = X(j);
    Y0 = Y(j);
    g = @(t) F(t,X0,Y0,pi/4); 
    q = quad(@(t)g(t),-2,2);
    g = @(t) F(t,X0,Y0,-5*pi/12); 
    q = q + quad(@(t)g(t),-2,2);
    g = @(t) F(t,X0,Y0,-13*pi/12); 
    q = q + quad(@(t)g(t),-2,2);
    if (mod(j,100) == 0)
        fprintf('%i operations finished.\n', j)
    end
    
    if q > 0.8
        Zint(ind) = q;
        Xint(ind) = X0;
        Yint(ind) = Y0;
        ind = ind + 1;
    end
end
Xint = Xint(1,1:ind-1);
Yint = Yint(1,1:ind-1);
Zint = Zint(1,1:ind-1);

min(Zint)
max(Xint)

plot3(Xint,Yint,Zint,'b*')
grid on;
hold on;































