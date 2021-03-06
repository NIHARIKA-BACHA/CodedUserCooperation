clear all;
close all;

Ni = 40;
Np = 40;
d = 7;
w = 1:80;
SNR = -5:10;
pbar=0;
Pe = 0; %BSC crossover probability
error_in = 0:Ni; % Range of errors in input Ni bits
error_out = 0:Np; % Range of errors in output Np bits
Pb=0;

for x=1:length(SNR)
% EsN0 calculation
Ec1 = 1;
Ec2 = 1;
EbN0 = 10^(SNR(x)/10);
EsN0 = EbN0;
N0 = (1/EbN0)*(Ec1); % Noise variance

% Compute the Pe for first phase BSC channel error
%pbar=1/2*(1-sqrt((Ec1*1/N0)/(1+(Ec1*1/N0))));
%pbar=1/2*erfc(sqrt(Ec1/N0));
pbar=0.01;
Pe=(1-(1-2*pbar)^d)/2;



Pb_km=0;
Pb_km_avg=0;
for k=error_in(2) : error_in(Ni+1)
    Pb_m=0;
    Pb_avg=0;
    for m=error_out(1) : error_out(Np+1)
        N_k_m=0;
        % Computing the number of codewords with input weight k and output
        % weight m, defined by N_k_m
        p_j_k=0;
        for j = max(0,k+d-Ni): min(d,k)
            if mod(j,2) == 1
                p_j_k = p_j_k + (nchoosek(d,j) * nchoosek((Ni-d),(k-j)))/nchoosek(Ni,k);
            else p_j_k = p_j_k;
            end
        end
        %p=sum(p_j_k)
        N_k_m = nchoosek(Ni,k)*nchoosek(Np,m)* p_j_k^m *(1-p_j_k)^(Np-m);
        
        % Computing the probability of m errors in Np parity check bits at
        % the relay node with e errors in Ni input bits
        P_m=0;
        N_kmtilda=0;
        
        for ktilda = error_in(2) : error_in(Ni+1)
            sum_d_kmtilda=0;
            N_km_tilda=0;
            for mtilda = error_out(1) : error_out(Np+1)
%                 if ktilda+mtilda <= Np
%                     d_kmtilda = randi(ktilda+mtilda);
%                 elseif ktilda+mtilda > Np
%                     d_kmtilda = randi(Np);
%                 end
%             sum_d_kmtilda = sum_d_kmtilda + d_kmtilda;
                p_jk=0;
                for j = max(0,ktilda+d-Ni): min(d,ktilda)
                    if mod(j,2) == 1
                         p_jk = p_jk + (nchoosek(d,j) * nchoosek((Ni-d),(ktilda-j)))/nchoosek(Ni,ktilda);
                    else p_jk = p_jk;
                    end
                 end
                N_km_tilda = N_km_tilda + nchoosek(Ni,ktilda)*nchoosek(Np,mtilda)* p_jk^mtilda *(1-p_jk)^(Np-mtilda);
            end
                        
%             if ktilda+m <= Np
%                 d_ktildam = randi(ktilda+m);
%             elseif ktilda+m > Np
%                 d_ktildam = randi(Np);
%             end
                pjk=0;
                for j = max(0,ktilda+d-Ni): min(d,ktilda)
                    if mod(j,2) == 1
                         pjk = pjk + (nchoosek(d,j) * nchoosek((Ni-d),(ktilda-j)))/nchoosek(Ni,ktilda);
                    else pjk = pjk;
                    end
                 end
            N_kmtilda = nchoosek(Ni,ktilda)*nchoosek(Np,m)* pjk^m *(1-pjk)^(Np-m);
            P_m = P_m + (N_kmtilda/N_km_tilda) * nchoosek(Ni,ktilda) * Pe^ktilda * (1-Pe)^(Ni-ktilda);
        end
        %P_mtotal = sum(P_m)
        
        Q_func = 1/2* erfc(((k*sqrt(Ec1)+m*sqrt(Ec2))/(sqrt(((k+m)*N0)/2)))/sqrt(2));
        
        Pb_m = Pb_m + N_k_m * (k/Ni) * P_m * Q_func;
        Pb_avg = Pb_avg + N_k_m * (k/Ni) * Q_func;
    end
    Pb_km = Pb_km + Pb_m;
    Pb_km_avg = Pb_km_avg + Pb_avg;
end
Pb(x) = Pb_km
Pbavg(x) = Pb_km_avg
end
semilogy(SNR,Pb, 'b-*'), hold on,
semilogy(SNR,Pbavg, 'b-'),
legend('P_b with error','Avg. P_b without error'),
xlabel('SNR(dB)'),
ylabel('P_{b}'),
%axis([SNR(1) max(SNR) 1e-7 1e0]);
