A = load('probs.txt');
x = 1: 1: size(A, 1);
plot(x, A(:, 2), 'b.');
grid on;