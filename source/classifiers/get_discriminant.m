function h = get_discriminant(X)
% get_discriminaton returns the the
% value of discriminant function
% for a vector X. Discriminant function
% is being calculated with fixed parameters
% for both classes, i.e. M1,S1,M2,S2

M1 = [4;8];
S1 = [2 -0.35; -0.35 1.9];
M2 = [4;4];
S2 = [0.9 0.7; 0.7 0.9];

M3 =[-4;6];
S3 = [1.5 0.5; 0.5 1.5];

f11 = 1/(2*pi)^(2/2)/det(S1)^0.5*exp(-0.5*(X-M1)'*S1^(-1)*(X-M1));
f12 = 1/(2*pi)^(2/2)/det(S2)^0.5*exp(-0.5*(X-M2)'*S2^(-1)*(X-M2));
f1 = 0.6*f11+0.4*f12;
f2 = 1/(2*pi)^(2/2)/det(S3)^0.5*exp(-0.5*(X-M3)'*S3^(-1)*(X-M3));
h = log(f2) - log(f1);
end

