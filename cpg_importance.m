b = h5read('weights.h5', '/dense_3/dense_3/bias:0') + h5read('weights.h5', '/dense_3/dense_3/kernel:0')*h5read('weights.h5', '/dense_2/dense_2/bias:0') + h5read('weights.h5', '/dense_3/dense_3/kernel:0')*h5read('weights.h5', '/dense_2/dense_2/kernel:0')*h5read('weights.h5', '/dense_1/dense_1/bias:0');
w = h5read('weights.h5', '/dense_3/dense_3/kernel:0') * h5read('weights.h5', '/dense_2/dense_2/kernel:0') * h5read('weights.h5', '/dense_1/dense_1/kernel:0');
for i = 1:18
    [k(i,:), p(i,:)] = sort(w(i,:), 'descend');
end
kk = k(:, 1:1000);
pp = p(:, 1:1000);
csvwrite('cpg_importance.csv',pp)