function [br1,br2,eps_vector] = wald_sequential(N,knum)
% wald_sequential implements Wald's sequential
% test for classes 1 i 2, based on vectors X1 and X2.
% During test mean number of samples necessary for
% making a decision is being calculated for fixed
% values of Eps1 or Eps2 with respect to the 
% progressive changes of the other parameter.
% The function returns vectors which represent
% average number of samples needed for both
% classes, with used error vector as well.

M1 = [4;8];
S1 = [2 -0.35; -0.35 1.9];
M2 = [4;4];
S2 = [0.9 0.7; 0.7 0.9];

M3 =[-4;6];
S3 = [1.5 0.5; 0.5 1.5];

f11 = []; f12 = []; 
f1 = []; f2 = []; 
hW = zeros(1,N);    % discriminant function
last_index = 0;     % index of last classified
br1 = [];   % number of necessary samples for class 1
br2 = [];   % number of necessary samples for class 1
eps_vector = logspace(-10,-2,100);

Eps2 = 1e-5;
for eps = eps_vector
 
    hW = [];
    Eps1 = eps; 
    A = (1-Eps1)/Eps2; a = -log(A); 
    B = Eps1/(1-Eps2); b = -log(B);
    
    broj1 = [];
    last_index = 0;
    [X1,X2] = class_generator(N);
    
    if knum==1
        Xchosen = X1;
    else
        Xchosen = X2;
    end
    
    for i=1:N
        X = Xchosen(:,i);
        hW(i) = get_discriminant(X);
        Sm = sum(hW(last_index+1:i));
        SS(i) = Sm;
        
        if Sm >= b
            broj1(end+1) = i-last_index;
            last_index = i;
        elseif Sm <= a
            broj1(end+1) = i-last_index;
            last_index = i;
        else
            if i==N
                break
            end
        end
    end
    br1(end+1) = mean(broj1);
end

Eps1 = 1e-5;
for eps = eps_vector
    
    hW = [];
    Eps2 = eps;
    A = (1-Eps1)/Eps2; a = -log(A); 
    B = Eps1/(1-Eps2); b = -log(B);
    
    broj2 = []; 
    last_index = 0;
    [X1,X2] = class_generator(N);
    
    if knum==1
        Xchosen = X1;
    else
        Xchosen = X2;
    end
    
    for i=1:N
        X = Xchosen(:,i);
        hW(i) = get_discriminant(X);
        Sm = sum(hW(last_index+1:i));

        if Sm >= b
            broj2(end+1) = i-last_index;
            last_index = i;
        elseif Sm <= a
            broj2(end+1) = i-last_index;
            last_index = i;
        else
            if i==N
                break
            end
        end
    end
    br2(end+1) = mean(broj2);
end

end

