clear all, close all, clc
%% Generating classes 

N = 500; % number of samples per class

% Generating class A samples
TetaA = rand(1,N)*2*pi;
MA = [0;0];
RA = 1.5;
Ra = RA*rand(1,N);
A0 = [Ra.*cos(TetaA); Ra.*sin(TetaA)] + MA*ones(1,N);

% Generating class B samples
TetaB = rand(1,N)*2*pi;
MB = [0;5];
RB = 2;
Rb = RB*rand(1,N);
B0 = [Rb.*cos(TetaB); Rb.*sin(TetaB)] + MB*ones(1,N);

% Generating class C samples
TetaC = rand(1,N)*2*pi;
MC = [4;0];
RC = 1.3;
Rc = RC*rand(1,N);
C0 = [Rc.*cos(TetaC); Rc.*sin(TetaC)] + MC*ones(1,N);

% Generating class D samples
TetaD = rand(1,N)*2*pi;
MD = [4;7];
RD = 2;
Rd = RD*rand(1,N);
D0 = [Rd.*cos(TetaD); Rd.*sin(TetaD)] + MD*ones(1,N);

   figure
plot(A0(1,:),A0(2,:),'rx',B0(1,:),B0(2,:),'bo',...
     C0(1,:),C0(2,:),'g*',D0(1,:),D0(2,:),'k.')
title('Class data')
legend('Class A', 'Class B', 'Class C','Class D','Location','East')
axis equal
hold off


%% K-means method implemented using elbow method

Z = [A0 B0 C0 D0];
num = 10;
J = zeros(1,num);

for K = 1:num
    
    clusterization = randi(K,2,4*N);
    clusterization(1,:) = realmax*ones(1,4*N);
    
    if K == 4
        A =  Z(:,find(clusterization(2,:)==1));
        B =  Z(:,find(clusterization(2,:)==2));
        C =  Z(:,find(clusterization(2,:)==3));
        D =  Z(:,find(clusterization(2,:)==4));  
        
        figure
        plot(A(1,:),A(2,:),'rx',B(1,:),B(2,:),'bo',...
             C(1,:),C(2,:),'g*',D(1,:),D(2,:),'k.')
        title('Initial clusterization L=4')
        legend('Class A', 'Class B', 'Class C','Class D','Location','East')
        axis equal
        hold off
    end
    
    for i=1:K
        elements = Z(:,find(clusterization(2,:)==i));
        cluster_centers(:,i) = mean(elements,2);
    end
    
    iter_max = 100; % max number of iterations
    iter = 1;       % counter of iterations
    done = 0;       % end of clusterization flag

    while (iter<=iter_max) && (done == 0)
        done = 1;
        for i=1:4*N
            sample = Z(:,i);
            
            for j=1:K
                distances(j) = ((sample-cluster_centers(:,j))'*...
                        (sample-cluster_centers(:,j))).^0.5;
            end
            [dmin,ind] = min(distances);
            
            if ind~=clusterization(2,i) 
                clusterization(:,i) = [dmin; ind];
                done = 0;
            else
                clusterization(1,i) = dmin;
            end
            
            for j=1:K
                if ~ismember(j,clusterization(2,:))
                    pos = (j-1)*5+1;
                    clusterization(2,pos:pos+4) = [j,j,j,j,j];
                end
            end
        end
        
        % Clusters centers
        for j=1:K
            elements = Z(:,find(clusterization(2,:)==j));
            cluster_centers(:,j) = mean(elements,2);
        end
        iter = iter + 1;
    end
   
    if K==4
        A =  Z(:,find(clusterization(2,:)==1));
        B =  Z(:,find(clusterization(2,:)==2));
        C =  Z(:,find(clusterization(2,:)==3));
        D =  Z(:,find(clusterization(2,:)==4));
        
        figure
        plot(A(1,:),A(2,:),'rx',B(1,:),B(2,:),'bo',...
             C(1,:),C(2,:),'g*',D(1,:),D(2,:),'k.')
        title('K-means clusterization result L=4')
        legend('Class A', 'Class B', 'Class C','Class D','Location','East')
        axis equal
        hold off
        
    elseif K==3 || K==5
        figure, hold all
        elements = [];
        for k=1:K
            elements =  Z(:,find(clusterization(2,:)==k));
            scatter(elements(1,:),elements(2,:))
            hold on
        end
        hold off
        axis equal
        title(['Number of clusters L=',num2str(K)])
    end

J(K) = mean(clusterization(1,:));
end

figure,
plot(1:num,J,'r-o','MarkerFaceColor','r')
grid on
title('Elbow method')
xlabel('Number of clusters')
ylabel('J')