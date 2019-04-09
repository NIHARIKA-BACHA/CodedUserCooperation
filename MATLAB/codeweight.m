


clear all;
close all;

Ni = 200; %Length of information bits vector
Np = 200; %Length of parity check bits vector
d = 9; % No. of 1's in a column of G matrix
N=[];
w = 1:Ni+Np;

%input = randn(1,Ni)<0.5;
%k= length(find(input));
%k=10;
Pjk=[];  
for k=1:Ni
    p=0;
    for j = max(0,k+d-Ni): min(d,k)
    if mod(j,2) ~= 0
        p = p + ((nchoosek1(d,j) * nchoosek1((Ni-d),(k-j)))/nchoosek1(Ni,k));
    end
    end
    Pjk(k)=p; %Prof. of j=1, given k 1's in the information vector
end
%Pjk
%p=sum(p_j_k);

for k=1:Ni
    for m=1:Np
     N(k,m) = nchoosek1(Ni,k)*nchoosek1(Np,m)*Pjk(k)^m*(1-Pjk(k))^(Np-m); %Total no. of codewords with input (information bits) weight k, and output (parity bits) weight m 
    end
end
%N
%figure(1)
%mesh(N)
N_w=zeros(1,length(w));
for i=1:length(w)
    N_sum=0;
 
  for k=1:Ni   
    for m=1:Np
        if k+m == i
            N_sum = N_sum + N(k,m);
        end
    end
  end
  N_w(i)=N_sum; % Total no. of codewords with weight, w
end
%N_w
figure(2)
plot(w,log(N_w)),  grid on,
title(['Code Weight Distribution of the Coded Cooperation Scheme, Ni=', num2str(Ni), ', Np=', num2str(Np),', d=',num2str(d)]),
xlabel('Weight, w'),
ylabel('No. of codewords with weight w, N(w)'),
%axis([0 length(w) min(N_w)-10 max(N_w)+10]);

