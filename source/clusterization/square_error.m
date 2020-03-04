clear all, close all, clc
%% Generating classes

N = 500;    % num of samples per class

% Generating class A samples
TetaA = rand(1,N)*2*pi;
MA = [5;5];
RA = 1.32;
Ra = RA*rand(1,N);
A0 = [Ra.*cos(TetaA); Ra.*sin(TetaA)] + MA*ones(1,N);

% Generating class B samples
TetaB = rand(1,N)*2*pi;
MB = [5;5];
RB = 2.75;
d = 0.9;
Rb = RB + d*rand(1,N);
B0 = [Rb.*cos(TetaB); Rb.*sin(TetaB)] + MB*ones(1,N);

% Class visuelization
figure, hold all
scatter(A0(1,:),A0(2,:),'rx')
scatter(B0(1,:),B0(2,:),'bo')
hold off
axis equal
legend('Class A', 'Class B','Location','NorthEast')
title('Classes')

%% Initial clusterization

Z = [A0 B0];
z = rand(2*N,1);
A = []; B = [];

Mpod = mean(Z,2);
for i=1:2*N
    if z(i)<0.5
         A = [A Z(:,i)];
     else
         B = [B Z(:,i)];
     end
end

% calculating initial clusters centers
MA = mean(A,2);
MB = mean(B,2);

SA = cov(A');
SB = cov(B');

% Initial clusterization visuelization
figure
plot(A(1,:),A(2,:),'rx',B(1,:),B(2,:),'bo')
axis equal
legend('Class A', 'Class B','Location','NorthEast')
title('Initial clusterization')
hold off

%% Square error method

iter_max = 100;
iter = 1;
done = 0;

while (iter<=iter_max) && (done==0)
    A1 = []; B1 = [];
    done = 1;
    
    for i=1:length(A)
        PA = length(A)/(2*N);
        PB = length(B)/(2*N);
        
        % calculating possible change of crit. function J
        % in case of possible reclusterization
        JA = 1/2*(A(:,i)-MA)'*(inv(SA))*(A(:,i)-MA) ...
             + 1/2*log(det(SA)) - 1/2*log(PA);
             
        JB = 1/2*(A(:,i)-MB)'*(inv(SB))*(A(:,i)-MB) ...
             + 1/2*log(det(SB)) - 1/2*log(PB);
         
         Jmin = min([JA, JB]);
         
         if Jmin==JA
             A1 = [A1 A(:,i)];
         else
             B1 = [B1 A(:,i)];
             done = 0;
         end
    end
    for i=1:length(B)
        PA = length(A)/(2*N);
        PB = length(B)/(2*N);
        
        % calculating possible change of crit. function J
        % in case of possible reclusterizatione
        JA = 1/2*(B(:,i)-MA)'*(inv(SA))*(B(:,i)-MA) ...
             + 1/2*log(det(SA)) - 1/2*log(PA);
             
        JB = 1/2*(B(:,i)-MB)'*(inv(SB))*(B(:,i)-MB) ...
             + 1/2*log(det(SB)) - 1/2*log(PB);
         
         Jmin = min([JA, JB]);
         
         if Jmin==JA
             A1 = [A1 B(:,i)];
             done = 0;
         else
             B1 = [B1 B(:,i)];
         end
    end
    
    % preventing index error in case that
    % algorith fails to detect some classes
    if isempty(A1)
        A1 = B1(:,1:2);
        B1 = B1(:,3:end);
    end
    
    if isempty(B1)
        B1 = A1(:,1:2);  
        A1 = A1(:,3:end);
    end
    
    clear A B
    A = A1; B = B1;
    clear A1 B1
    
    MA = mean(A,2); SA = cov(A');
    MB = mean(B,2); SB = cov(B');
    iter = iter + 1;
end

MA = mean(A,2);
MB = mean(B,2);

figure
plot(A(1,:),A(2,:),'rx',B(1,:),B(2,:),'bo')
title('Square decomposition clusterization result')
legend('Class A', 'Class B','Location','NorthEast')
axis equal
hold off