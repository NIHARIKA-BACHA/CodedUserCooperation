function Pe = cuc_errorbound2() %#codegen

clear all;
close all;

Ni = 50;
Np = 50;
d = 7;
w = 1:Ni+Np;
SNR = 0:1:15;
pbar=0;
pe = 0; %BSC crossover probability
%error_in = 0:Ni; % Range of errors in input Ni bits
%error_out = 0:Np; % Range of errors in output Np bits
Pe=zeros(1,length(SNR));
Pe1=0;

%Codeweight distribution
N = codeweight(Ni, Np, d);
%N
for x=1:length(SNR)
% EsN0 calculation
Es1 = 1;
Es2 = 1;
R_LDGM = Np/(Ni+Np);
Ec1 = Es1;
Ec2 = Es2;
EbN0 = 10^(SNR(x)/10);
EsN0 = EbN0*R_LDGM;
N0coded = (1/EbN0)*(Ec2/R_LDGM); % Noise variance LDGM
N0uncoded = (1/EbN0)*Ec1; % Noise variance

% Compute the Pe for first phase BSC channel error
%pbar=1/2*(1-sqrt((Ec1*1/N0uncoded)/(1+(Ec1*1/N0uncoded))));
pbar= 1/2* erfc((sqrt(2*Ec1/N0uncoded))/sqrt(2));
%pbar=1/2*(1-sqrt(EbN0/(1+EbN0)));
%pe=1/2*erfc(sqrt(Ec1/N0uncoded));
%pbar=0.0001;
pe=(1-(1-2*pbar)^d)/2;
%pe = 0.0001;

%Pairwise error probability
P = pairwise(Ni, Np, N0coded, Ec1, Ec2, pe, N);


Pe1=0;
for k=1:Ni
    for m=1:Np
      %Pe1
      sum_nkm=0;
      for nkm=1:Np
          sum_nkm = sum_nkm + N(k,nkm);
      end
      Pe1 = Pe1 + P(k,m)*(1/2* erfc(sqrt((2*k*Ec1)/(2*N0uncoded))));
      %Pe1 = Pe1 + ((1/(Ni+Np))*P(k,m));
    end
   
end
Pe(x) = Pe1;
fprintf('\nSimulation done for SNR: %d\n',SNR(x));
end

% % Result storage
% save_mat_name = ['./results/UnionBound_Ni' int2str(Ni) 'Np' int2str(Np) 'd' int2str(d)];
% save (save_mat_name,'Pe');
% % Loading saved results for plotting
% load_mat_name = ['./results/UnionBound_Ni' int2str(Ni) 'Np' int2str(Np) 'd' int2str(d)];
% load(load_mat_name);
% Pe = Pe;
    
figure(4)
semilogy(SNR,Pe, 'b-*'), hold on, grid on;
%semilogy(SNR,Pbavg, 'b-');
%legend('P_e with error','Avg. P_e without error');
xlabel('SNR(dB)');
ylabel('P_{e}');
%axis([SNR(1) max(SNR) 1e-7 1e0]);
end

% Function to compute pairwise error probability
    function P = pairwise(Ni, Np, N0, Ec1, Ec2, pe, N)
    
    count =0;
    count1=0;
    P=zeros(Ni,Np);
    Prob_Nsum = zeros(Ni,Np);
    for k=1:Ni
      for m=1:Np    
        PEP = 0;
        Prob_N = zeros(Ni,Np);
        Prob_Csum = zeros(Ni,Np);
        dE = 2*sqrt(k*Ec1 + m*Ec2);
                       
        for e=1:Np
            
            %Sign of b
            %theta = acos(k*sqrt(e*Ec1) / sqrt(k*Ec1 + m*Ec2));
            theta = acos(m*sqrt(e*Ec2) / sqrt(k*Ec1 + m*Ec2));
            %theta = acos(m*sqrt(Ec1) / sqrt(k*Ec1 + m*Ec2));
            %a = sqrt(e*(2*sqrt(Ec1))^2);
            

            if abs(theta) <= pi/2
                signb = 1;
            elseif abs(theta) > pi/2 && abs(theta) < 3*pi/2
                signb = -1;
            elseif abs(theta) >= 3*pi/2
                signb = 1;
            end
            if signb == -1, count = count + 1; else count1 = count1+1; end
            %signb
            
            b =  signb * 2*e*m*Ec1/sqrt(k*Ec1+m*Ec2);
            %b = signb * 2*e*m*sqrt(Ec1*Ec2)/sqrt(k*Ec1+ m*Ec2);
            %abs(theta)
            

            %dE/2-b
            %Q function for (dE/2 - b) as a reference point
            %Q_func= 1/2* erfc((sqrt(2/N0)*((2*(k*Ec1+m*Ec2) - (4*e*k*Ec1)*signb)/(2*sqrt(k*Ec1+m*Ec2))))/sqrt(2));
            %Q_func= 1/2* erfc((sqrt(2/N0)*(dE/2)-b)/sqrt(2));
            Q_func= 1/2* erfc((sqrt(2/N0)*(((k*Ec1+m*Ec2-2*e*m*sqrt(Ec1*Ec2))*signb)/sqrt(k*Ec1+m*Ec2)))/sqrt(2));
            
            %Q function for (dE/2) as a reference point
            %Q_func= 1/2* erfc((sqrt(2*(k*Ec1+m*Ec2)/N0))/sqrt(2));
            
            %Q function for (dE - b) as a reference point
            %Q_func= 1/2* erfc((sqrt(2/N0)*((2*(k*Ec1+m*Ec2) - e*m*sqrt(Ec1*Ec2))/(sqrt(k*Ec1+m*Ec2))))/sqrt(2));
            Prob_C = zeros(Np,Np);
                        
            sum_mddash=0;
            for mddash=1:Np
                
                sum_mddash = sum_mddash + N(e, mddash);
            end
            %sum_mddash
            
            sum_mdash = 0;
            for mdash=1:Np
            
                sum_ii=0;
                for ii=max(0,m+mdash-Np):min(m,mdash)%0:mdash 
%                     m
%                     ii
%                     Np-m
%                     mdash-ii
                    %if max(mdash,m) >= ii && Np-max(mdash,m) >= mdash-ii && mdash > ii
                    %if ii <= m && mdash-ii <= Np-m
                        sum_ii = sum_ii + (nchoosek1(m, ii) * nchoosek1(Np-m, mdash-ii));
                    %end
                end
                %sum_ii
                

                sum_i=0;
                for i=max(0,m+mdash-Np):min(m,mdash)%1:mdash 
                
                    %if max(mdash,m) >= i && Np-max(mdash,m) >= mdash-i && mdash > i
                    %if i <= m && mdash-i <= Np-m
                        sum_i = sum_i + (((nchoosek1(m, i) * nchoosek1(Np-m, mdash-i)) / sum_ii) * Q_func);
                        %sum_i = sum_i + ( Q_func );
                        %Prob_C(mdash,i) = (nchoosek1(m, i) * nchoosek1(Np-m, mdash-i)) / sum_ii;
                    %else Prob_C(mdash,i) = 0;
                    %end
                end
                
                %Prob_C
                %sum_i
                %sum_mdash = sum_mdash + sum_i;
                sum_mdash = sum_mdash + (sum_i * (N(e, mdash)/sum_mddash)); %Include the N/sum_N term or not               
                Prob_N(m,mdash) = N(m, mdash)/sum_mddash;
            end
            %sum_mdash
            %Q_func
            %Prob_Csum(e,:) = sum(Prob_C,1)./Np;
%             figure(1)
%             semilogy(1:length(Prob_Csum),Prob_Csum, 'g-*'), hold on, grid on,
%             figure(2)
%             semilogy(1:Np,Prob_N, 'r-*'), hold on, grid on,
            %Prob_N(e)=sum(Prob_N)/mdash;
            %sum_mdash
            PEP = PEP + (nchoosek1(Ni, e)* pe^e * (1-pe)^(Ni-e) *  sum_mdash);
        end
        %PEP
        P(k,m) = PEP; % Compute PEP for each k and m and store in P.
        %P(k,m) = 1/2* erfc((k*sqrt(Ec1)+m*sqrt(Ec1))/sqrt(((k+m)*N0)/2)/sqrt(2)); % Average BER for AWGN without inter-user channel error
        %Prob_Csum1(k,:) = sum(Prob_Csum,1)./Ni;
        %Prob_Nsum(k,:) = sum(Prob_N,1)./Ni;
%         figure(1)
%         semilogy(1:length(Prob_Csum1),Prob_Csum1, 'g-*'), hold on, grid on,
%         figure(2)
%         semilogy(1:Np,Prob_Nsum1, 'r-*'), hold on, grid on,
      end
    end
    %P
    %Prob_Csum2 = sum(Prob_Csum,1)./Ni;
    %Prob_Nsum2 = sum(Prob_Nsum,1)./Ni;
    %figure(1)
    %semilogy(1:length(Prob_Csum),Prob_Csum2, 'g-*'), hold on, grid on,
    %figure(2)
    %semilogy(1:length(Prob_Nsum),Prob_Nsum2, 'r-*'), hold on, grid on;
    %count
    %count1
    
%     w = 1:Ni+Np;
%     P_w=zeros(1,length(w));
%     for i=1:length(w)
%         P_sum=0;
% 
%       for k=1:Ni   
%         for m=1:Np
%             if k+m == i
%                 P_sum = P_sum + P(k,m);
%             end
%         end
%       end
%       P_w(i)=P_sum; % Total no. of codewords with weight, w
%     end
%     %N_w
%     figure(3)
%     semilogy(w,P_w),  grid on;
%     %title(['Code Weight Distribution of the Coded Cooperation Scheme, Ni=', num2str(Ni), ', Np=', num2str(Np),', d=',num2str(d)]),
%     xlabel('Weight, w');
%     ylabel('No. of codewords with weight w, N(w)');
    end

% Function to compute the distribution of codewords
    function N = codeweight(Ni, Np, d)

        N=zeros(Ni,Np);
        Pjk=zeros(1,Ni);  

        for k=1:Ni
            p=0;
            for j = max(0,k+d-Ni): min(d,k)
                if mod(j,2) ~= 0
                    p = p + ((nchoosek1(d,j) * nchoosek1((Ni-d),(k-j)))/nchoosek1(Ni,k));
                end
            end
            Pjk(k)=p; %Prob. of j=1, given k 1's in the information vector
        end
        %Pjk
        for k=1:Ni
            for m=1:Np
                N(k,m) = nchoosek1(Ni,k)*nchoosek1(Np,m)*Pjk(k)^m*(1-Pjk(k))^(Np-m); %Total no. of codewords with input (information bits) weight k, and output (parity bits) weight m 
            end
        end
    end


