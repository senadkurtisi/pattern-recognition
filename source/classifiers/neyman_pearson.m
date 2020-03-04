function [ T_NP, h_opt, Eps2 ] = neyman_pearson(f1, f2, Eps1, Eps2,...
                                               x_len, y_len)
% NEYMAN_PEARSON determines the threshold of
% the Neyman-Pearson classifier=-ln(mi) with
% respect to desired CFAR: Eps0 = Eps + Eps2

mi_vec = 0:0.01:5;
h_opt = [];

format long
Eps0 = Eps1+Eps2 % CFAR
Eps2 = 0;

for mi = mi_vec
    Eps2 = 0;
    for i = 1:x_len
        for j = 1:y_len
            hNP(i,j) = log(f2(i,j)) - log(f1(i,j))-log(mi);
            if hNP(i,j) < 0
                Eps2 = Eps2 + f2(i,j)*0.1*0.1; 
            end
            
        end
    end
    
    if Eps2>=Eps0
        mi_opt = mi;
        h_opt = hNP;
        break;
    end
end

T_NP = -log(mi_opt);
fprintf('Achieved type 2 error: %.8f\n', Eps2)

end

