function [X1, X2] = class_generator(N)
% class_generator is generating two classes
% with respect to the argument N which is
% the number of samples needed. Classes are
% being generated with fixed distribution
% parameters.

M1 = [4;8];
S1 = [2 -0.35; -0.35 1.9];
M2 = [4;4];
S2 = [0.9 0.7; 0.7 0.9];

M3 =[-4;6];
S3 = [1.5 0.5; 0.5 1.5];

% T = F*L^0.5 
[Fi1, L1] = eig(S1); T1 = Fi1*L1^0.5;
[Fi2, L2] = eig(S2); T2 = Fi2*L2^0.5;
[Fi3, L3] = eig(S3); T3 = Fi3*L3^0.5;

%K1
X1 = zeros(2,N);
p = rand(1,N);
pnom = randn(2,N);
X1 = (p<0.6).*(T1*pnom+M1*ones(1,N))+...
 (p>0.6).*(T2*pnom+M2*ones(1,N));

%K2
X2 = zeros(2,N);
X2 = T3*randn(2,N)+M3*ones(1,N);

end
