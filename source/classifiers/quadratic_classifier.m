clear all, close all, clc
%% Quadratic classifier

N = 1000;
% Class 1: full circle
% Class 2: Torus/Ring

M1 = [8;8];
R1 = 2;
TetaX = rand(1,N)*2*pi;
Rx = R1*rand(1,N);
X = [Rx.*cos(TetaX); Rx.*sin(TetaX)] + M1*ones(1,N);

M2 = [8;8];
R2 = 3; d=2;
TetaY = rand(1,N)*2*pi;
Ry = R2+d*rand(1,N);
Y = [Ry.*cos(TetaY); Ry.*sin(TetaY)] + M2*ones(1,N);

figure
plot(X(1,:),X(2,:),'r.', Y(1,:), Y(2,:),'bo')
title('Non-linearly separable classes')
legend('Class 1', 'Class 2')
axis equal

%% Quadratic classifier creation
Gama = ones(2*N,1);
Gama(N+1:end) = ones(N,1)*1.5;

Z = [ ones(1,N) -ones(1,N); ...
    X -Y;...
    X(1,:).^2 -Y(1,:).^2; ...
    2*X(1,:).*X(2,:) -2*Y(1,:).*Y(2,:);...
    X(2,:).^2 -Y(2,:).^2
    ];
    
W = pinv(Z')*Gama;
v0 = W(1); v = W(2:3); Q = [W(4:5); W(5:6)];

syms xp yp
[xp, yp, ~,~] = solve(v0+v(1)*xp+v(2)*yp+...
    xp^2*Q(1)+xp*yp*2*Q(2)+yp^2*Q(4), xp, yp, 'returnconditions',true);

z = 4:0.0001:12;
xp = eval(xp);  xplot = [xp(1,:) fliplr(xp(2,:))];
yp = eval(yp);  yplot = [yp(1,:) fliplr(yp(2,:))];

xplot1 = xplot(imag(xplot)==0);
yplot1 = yplot(imag(xplot)==0);

figure(1)
hold on
plot(xplot1,yplot1,'k')
legend('X2','X1','Quadratic classifier')
