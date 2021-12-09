function T = Tensor_2x2bin(T)


A = T(1:2:end,1:2:end,:);
B = T(1:2:end,2:2:end,:);
C = T(2:2:end,1:2:end,:);
D = T(2:2:end,2:2:end,:);
T = (A+B+C+D)/4;
