clear all, close all, clc
%% Generating classes

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

N = 500; %num of elements per class

%K1
X1 = zeros(2,N);
p = rand(1,N);
pnom = randn(2,N);
X1 = (p<0.6).*(T1*pnom+M1*ones(1,N))+...
     (p>0.6).*(T2*pnom+M2*ones(1,N));

%K2
X2 = zeros(2,N);
X2 = T3*randn(2,N)+M3*ones(1,N);

% Class visualization
figure
plot(X2(1,:),X2(2,:),'rx')
axis equal
hold on
plot(X1(1,:),X1(2,:),'bo')
hold off
legend('Class 2','Class 1')
title('Classes')

%% 
x = -8:0.1:8;
y = -1:0.1:11;

for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i);y(j)];
        f11(i,j) = 1/(2*pi)^(2/2)/det(S1)^0.5*...
                    exp(-0.5*(X-M1)'*S1^(-1)*(X-M1));
        f12(i,j) = 1/(2*pi)^(2/2)/det(S2)^0.5*...
                    exp(-0.5*(X-M2)'*S2^(-1)*(X-M2));
        f1(i,j) = 0.6*f11(i,j)+0.4*f12(i,j);
        f2(i,j) = 1/(2*pi)^(2/2)/det(S3)^0.5*...
                    exp(-0.5*(X-M3)'*S3^(-1)*(X-M3));
        
        h(i,j) = log(f2(i,j))-log(f1(i,j));
    end
end

%% d2 

dk = [1.2 4 10];
for d2 = dk
    teta = 0:0.01:2*pi;
    Xpom = [d2^0.5*cos(teta); d2^0.5*sin(teta)];
    Y2=T3*Xpom + M3*ones(1,length(teta));
    figure(1) 
    hold on
    plot(Y2(1,:),Y2(2,:),'k.')
end

%K2
f2max = max(f2(:));
prag1 = f2max*exp(-dk(1)/2);
prag2 = f2max*exp(-dk(2)/2);
prag3 = f2max*exp(-dk(3)/2);
contour(x,y,f2',[prag1 prag1],'m');
contour(x,y,f2',[prag2 prag2],'m');
contour(x,y,f2',[prag3 prag3],'m');

%K1
f1max = max(f1(:));
prag1 = f1max*exp(-dk(1)/2);
prag2 = f1max*exp(-dk(2)/2);
prag3 = f1max*exp(-dk(3)/2);

Y1 = [];
for i=1:length(x)
    for j =1:length(y)
        if abs(f1(i,j)-prag1)<5e-3
        Y1 = [Y1 [x(i);y(j)]];
        end
    end
end

figure(1), hold on
plot(Y1(1,:),Y1(2,:),'k.')

Y1 = [];
for i=1:length(x)
    for j =1:length(y)
        if abs(f1(i,j)-prag2)<5e-4
        Y1 = [Y1 [x(i);y(j)]];
        end
    end
end
figure(1), hold on
plot(Y1(1,:),Y1(2,:),'m.')

Y1 = [];
for i=1:length(x)
    for j =1:length(y)
        if abs(f1(i,j)-prag3)<1e-4
        Y1 = [Y1 [x(i);y(j)]];
        end
    end
end
figure(1),hold on
plot(Y1(1,:),Y1(2,:),'g.')
legend('Class 2','Class 1')

%% Hypothesis tests

P2 = 0.5; P1 = 0.5;
C12 = 3;  C21 = 10;
C11 = 0;  C22 = 0;

Tb = log(P1/P2);    % bayes min. error probability treshold
Tc = log((C21-C11)*P1/((C12-C22)*P2)); % bayes min. price threshold

x = -8:0.1:8;
y = -1:0.1:11;
for i = 1:length(x)
    for j = 1:length(y)
        hB(i,j) = log(f2(i,j))-log(f1(i,j))-Tb;
        hC(i,j) = log(f2(i,j))-log(f1(i,j))-Tc;
    end
end

% Error estimation for Bayes min. error probability classificator
Eps1 = 0;   % type 1 error
Eps2 = 0;   % type 2 error

for i=1:length(x)
    for j=1:length(y)
        if hB(i,j) < 0 
            Eps2 = Eps2 + f2(i,j)*0.1*0.1; 
        else
            Eps1 = Eps1 + f1(i,j)*0.1*0.1;
        end
    end
end
fprintf("Type 1 error: %5.10f\n",Eps1)
fprintf("Type 2 error: %5.10f\n",Eps2)

%% Neyman-Pearson

[T_NP, hNP] = neyman_pearson(f1, f2, Eps1, Eps2, length(x), length(y)); 

figure(3), hold all
plot(X2(1,:),X2(2,:),'rx')
axis equal
plot(X1(1,:),X1(2,:),'bo')
contour(x,y,hNP',[T_NP T_NP],'k');
contour(x,y,hB',[Tb Tb],'m-.');
contour(x,y,hC',[Tc Tc],'g');
title('Two class classification')
legend('Class 2','Class 1','Neyman-Pearson','Min. error probability'...
       ,'Bayes min. price','Location','SouthWest')
txt = {['C_{11} = C_{22} = ', num2str(C11)], ...
       ['C_{12} = ', num2str(C12)], ...
       ['C_{21} = ', num2str(C21)]};
text(5, 1.65, txt, 'FontSize', 8.5, 'HorizontalAlignment','left')