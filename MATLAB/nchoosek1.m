

function Nk = nchoosek1(N, k) %#eml
%N = 200;
%k = 5;

i = N:-1:1;
ii = k:-1:1;
j = N-k:-1:1;

%Nk = exp(sum(log(i))) / (exp(sum(log(ii))) * exp(sum(log(j))));
Nk = exp(sum(log(i)) - (sum(log(ii)) + sum(log(j))));
end


