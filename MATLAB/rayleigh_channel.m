% time varying Rayleigh flat fading channel with Clark's model 
% using filter method
% Ns number of samples
function out=rayleigh_channel(Ns, N0)
%% parameters
%Ns=32;
%N0=0.1;
fc=2E9;
% velocity
v=50*1000/3600;
% speed of light
c=3E8;
% doppler frequency
fD=(fc/c)*v;
Ds=2*fD;%Doppler spread
% sampling time
Ts=.1e-3;                     
Fs=1/Ts;
% filter window length
Ng=10;
% autocorrelation window length
w=400;
%% Complex Gaussian input
x=sqrt(N0/2)*(randn(1,Ns)+1i*randn(1,Ns));
%x=sqrt(N0/2)* ones(1,Ns);
%% compute g(t) and g_hat(t)
t=(-Ng:Ng)*Ts;
for k=1:2*Ng+1
    if(k==Ng+1)
        g(k)=((pi*fD)^(1/4))/gamma(5/4);
    else
        g(k)= besselj(1/4,2*pi*fD*abs(t(k)))/((abs(t(k)))^(1/4));
    end
end
%% compute K such that ... 
g_hat=g./sqrt(sum(g.^2));
% verify K
sum_g_hat=sum(g_hat.^2);
%% channel impulse response
cn=conv(x,g_hat)*Ts;
%stem(cn); hold on
%% reject transient
cn=cn(2*Ng+1:length(cn)-2*Ng);
out=cn/sqrt(var(cn));
stem(out);
%plot(abs(out.^2));hold off
end