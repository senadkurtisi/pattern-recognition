clear all, close all, clc
%% Wald sequential test

N = 500;

for knum = [1,2]
    [br1, br2, eps_vector] = wald_sequential(N,knum);
    n = length(eps_vector);
    
    if knum==1
        kls = 'first';
    else
        kls = 'second';
    end
    
    figure, semilogx(eps_vector,br1,'Linewidth',2)
    title(['Number of samples needed for ', ...
            kls ,' class ',...
            char(949),'_{2}','=10^{-05}'])
    xlabel([char(949), '_{1}[%]'])

    figure, semilogx(eps_vector,br2,'Linewidth',2)
    title(['Number of samples needed for ', ...
            kls, ' class ',...
            char(949),'_{1}','=10^{-05}'])
    xlabel([char(949), '_{2}[%]'])
end