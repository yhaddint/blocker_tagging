clear;clc
% cs parameters
N_cs=900;
M_cs=N_cs*0.15;
% TW parameters
N=31;M=30;N_stage=2;K=31;

C_STPE = (45)*log(45)*2*4500/100;
C_TW = N*M*N_stage*K/2+M*N_stage*log2(N)*K+K*M*N_stage+log2(M)*K*N_stage;
%C_TW = M*N_stage*log2(N)*K+K*M*N_stage+log2(M)*K*N_stage;


K_interest = 30;
k=1:K_interest;
C_OMP_ite = 2*N_cs*M_cs+2*M_cs*k+2*k.^2+2*M_cs*k;

C_OMP=zeros(1,K_interest);
C_OMP(1) = C_OMP_ite(1);

for ii=2:K_interest
    C_OMP(ii)=C_OMP_ite(ii)+C_OMP(ii-1);
end

%%
figure
plot_setting
semilogy(1:K_interest, ones(1,K_interest)*C_TW,'r');hold on
semilogy(1:K_interest, ones(1,K_interest)*C_STPE,'g');hold on
semilogy(1:K_interest,C_OMP,'b');hold on

xlim([1,30])
legend('C2I2','STSA','CS-OMP')
grid on
ylabel('Number of Computations')
xlabel('Number of Actual Interferers')