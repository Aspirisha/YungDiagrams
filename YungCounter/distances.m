A = load('distances.txt');
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
colorbar
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'metricHist2.jpg');

%strict dists
A = load('StrictDistances.txt');
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
colorbar
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, 'metricHist2.jpg');