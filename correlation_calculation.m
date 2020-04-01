Dataset = ...
[6.28	4.45
6.21	4.64
6.4     4.89
6.44	5.02
6.89	5.51
6.72	5.39
6.6     5.3
6.47	5.26
5.26	4.08
5.47	4.29
6.65	5.51
6.09	5.04
6.84	5.78
6.81	5.77
6.83	5.85
6.76	5.98
6.72	6.02
6.87	6.19
6.35	5.69
6.36	5.75
6.57	6.02
6.55	6.04
6.48	6
6.33	5.87
6.36	5.95
6.37	6.02
6.59	6.31
6.46	6.35
5.56	5.59
6.18	6.27];

U3 = Dataset(:,2);
max_U3 = max(abs(U3));
U3 = U3./max_U3;

H = Dataset(:,1);
max_H = max(abs(H));
H = H./max_H;

% Pearson r
sumU3 = sum(abs(U3))/size(U3,1);
reference = abs(U3)-sumU3;
varreference = sqrt(var(reference));

% Spearman's rho
[huy_s,huy_i] = sort(abs(U3));
R_huy = zeros(size(U3,1),size(U3,2));

for rank = 1:size(huy_i,1)
    R_huy(huy_i(rank,1),1) = rank;
end

%Kendall's tau
T_huy = zeros(size(U3,1),size(U3,1));
for ver = 1:size(T_huy,1)
    for hor = ver:size(T_huy,2)
        T_huy(ver,hor) = sign(abs(U3(ver,1))-abs(U3(hor,1)));
    end
end


%calculate r
sumH = sum(abs(H))/size(H,1);
compare = abs(H)-sumH;
varcompare = sqrt(var(compare));

cross_correlation = sum(reference.*compare)/varreference/varcompare/30;
H_crr_error = cross_correlation;

display(['Pearson test correlation value:',num2str(H_crr_error)])

[hyb_s,hyb_i] = sort(abs(H));
R_hyb = zeros(size(H,1),size(H,2));

for rank = 1:size(hyb_i,1)
    R_hyb(hyb_i(rank,1),1) = rank;
end

%calculate rho
H_rho_error = 1-(6*sum((R_huy-R_hyb).^2))/(size(H,1)*(size(H,1)^2-1));
display(['Spearman test correlation value:',num2str(H_rho_error)])


%calculate tau
discordance = 0;
concordance = 0;
T_hyb = zeros(size(U3,1),size(U3,1));
for ver = 1:size(T_hyb,1)
    for hor = ver:size(T_hyb,2)
        T_hyb(ver,hor) = sign(abs(H(ver,1))-abs(H(hor,1)));
        D = abs(T_huy(ver,hor)-T_hyb(ver,hor));
        C = abs(T_huy(ver,hor)+T_hyb(ver,hor));
        if D ==2
            discordance = discordance+1;
        elseif C ==2
            concordance = concordance+1;
        end
    end
end


H_tau_error = (concordance-discordance)/(size(H,1)*(size(H,1)-1))*2;
display(['Kendall test correlation value:',num2str(H_tau_error)])
