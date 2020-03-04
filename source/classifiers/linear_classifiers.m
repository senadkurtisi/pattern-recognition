clear all, close all, clc

%% Generating two 2D classes

N = 1000;
M2 = [0;10];
S2 = [0.9 0.7; 0.7 0.9];
M1 = [4.2;8];
S1 = [0.5 -0.1; -0.1 0.5];

[Fi1, L1] = eig(S1); T1 = Fi1*L1^0.5;
[Fi2, L2] = eig(S2); T2 = Fi2*L2^0.5;

X1 = T1*randn(2,N)+M1*ones(1,N);

X2 = T2*randn(2,N)+M2*ones(1,N);

figure
plot(X2(1,:),X2(2,:),'r.')
axis equal
hold on
plot(X1(1,:),X1(2,:),'bo')
hold off
title('Class data')
legend('Class 2','Class 1')

%% II iterative method
Nerr_opt = N; % initially equal to number of samples
s_opt = 0;
v0_opt = 0;
E_vec = [];

for s = 0:0.01:1
    V = (s*S1 + (1-s)*S2)^-1*(M2-M1);
    sigma1 = V'*S2*V;
    sigma2 = V'*S2*V;
    y_1 = V'*X1;
    y_2 = V'*X2;
    
    % calculating the next loop limits {v0}
    y_min = V'*X1(:,1); y_max = y_min;
    for i=1:N
        if(V'*X1(:,i)<y_min) 
            y_min = V'*X1(:,i);
        end
        if(V'*X2(:,i)<y_min) 
            y_min = V'*X2(:,i);
        end
        if(V'*X1(:,i)>y_max) 
            y_max = V'*X1(:,i);
        end
        if(V'*X2(:,i)>y_max) 
            y_max = V'*X2(:,i);
        end
    end
   
    N_errors = 2*N;
    v0_temp = 0;
    for v0 = -y_max:(-y_min)
        E = 0;
        
        E = E + sum(y_1 > (-v0)) + sum(y_2 <(-v0));
        if E<N_errors
            N_errors = E;
            v0_temp = v0;
        end
    end
    
    if N_errors<Nerr_opt
        Nerr_opt = N_errors;
        v0_opt = v0_temp;
        s_opt = s;
    end
    
    E_vec(end+1) = Nerr_opt;
end

Vopt = (s_opt*S1+(1-s_opt)*S2)^(-1)*(M2-M1);
x= -4:0.1:10;
y= -5.4:0.1:14;

for i=1:length(x)
    for j=1:length(y)
        X=[x(i); y(j)];
        h(i,j)=Vopt'*X+v0_opt;
    end
end

figure(1)
hold on
contour(x,y, h', [0 0], 'g');
hold off
legend('X2','X1','II iterative method classifier','Location','SouthWest')

% calculating prefered plot limits
x_min = min(min(X1(1,:),X2(1,:)))-0.3;
x_max = max(max(X1(1,:),X2(1,:)))+0.3;
y_min = min(min(X1(2,:),X2(2,:)))-0.3;
y_max = max(max(X1(2,:),X2(2,:)))+0.3;
axis([x_min x_max y_min y_max])

%% Desired output method

% defining desired output
Gama = ones(2*N,1);
Gama(1:N) = ones(N,1)*2;

Z = [ones(1,N) ,-ones(1,N); X1 ,-X2];
W = pinv(Z')*Gama;

% extracting necessary parameters
% for creating classification fcn
v0 = W(1); v = [W(2);W(3)];


xp = [0 5.5];
yp= -v(1)/v(2)*xp - v0/v(2);

figure(1)
hold on 
plot(xp,yp,'k-.')
legend('X2','X1','II iterative method classifier',...
       'Desired output classifier','Location','SouthWest')

