% 2016/09/05
% LMS based predistortion with memory depth 2

clear;%clc;

%% all parameters

SNR_dB = 30;
SNR = 10^(SNR_dB/10);

EZmodel = 1;
% rand('seed',1)
% M = 4;
% s = exp(1j*rand(M,1)*2*pi)/sqrt(M/2);
% s = [1+j;-1]/2;
beta1_0 = 1.0513+1j*0.0904;
beta3_0 = -0.0542-1j*0.2900;
beta1_1 = -0.0680-1j*0.0023;
beta3_1 = 0.2234+1j*0.2317;

% beta1_0 = 1.0513;
% beta3_0 = -0.0542;
% beta1_1 = -0.0680;
% beta3_1 = 0;
beta = [beta1_0;beta3_0;beta1_1;beta3_1];
% y = s;



% for mm=1:M
%     phi(:,mm) = [s;s(1)*abs(s).^2;s(2)*abs(s).^2;s(3)*abs(s).^2;s(4)*abs(s).^2;];
% %     phi(:,mm) = [s;s(1)*abs(s).^2;s(2)*abs(s).^2];
% end

% NLAP model
% phi2 = zeros(M,2);
% phi2(:,1) = x;
% for mm=1:M
%     phi2(mm,2) = x(mm)*abs(x(mm))^2+2*x(mm)*(sum(abs(x).^2)-abs(x(mm))^2);
% end

% analog beamforming and combining
% Nt = 128;
% AOD = linspace(-pi/5,pi/5,M).';
% Frf = zeros(Nt,M);
% for mm=1:M
%     Frf(:,mm) = exp(1j*pi*(0:Nt-1).'*sin(AOD(mm)));
% end
% 

%% LMS iteration
% stepsize_pool = [2e-3,1e-3,0.5e-3,0.5e-3,0.5e-3,0.6e-3,0.7e-3];
% ite_num_pool = [5e3,5e3,5e3,7.5e3,1e4,1e4,1e4];


M_pool = 2:20;
results = zeros(3,length(M_pool));
for beam_index=1:length(M_pool)
% ite_num = ite_num_pool(beam_index);
% MSE = zeros(ite_num,1);
% stepsize = stepsize_pool(beam_index);
M = M_pool(beam_index)
P = ones(M,1)/M;
% noise_pow = ones(M,1)/SNR;
%initialization of vector, matrix
% s = zeros(M,ite_num*2);
x = zeros(M,1);
% h=zeros(ite_num,1);
% h=zeros(M,M);
% w=zeros(2*(M+M^3),M);
% for mm=1:M
%     w(mm,mm) = 1;
% end

%     s = exp(1j*rand(M,ite_num)*2*pi)/sqrt(M/2);

% for mm=1:M
%     s(mm,:) = ones(1,2*ite_num)*sqrt(P(mm));
%     s_raw = ((randn(1,2*ite_num))+1j*(randn(1,2*ite_num)))/sqrt(2);
% %     s_raw = exp(1j*rand(1,ite_num)*2*pi)/sqrt(2);
% %     s_raw = ((rand(1,2*ite_num)*2-3)+1j*(randn(1,2*ite_num)*2-3))/sqrt(2);
% %         s_raw = rand(1,ite_num)*3;
%     PAPR = 2.5;
%     index = find(abs(s_raw).^2<PAPR);
%     s(mm,index) = s_raw(index)*sqrt(P(mm));
% end
% 
% x = zeros(M,2*ite_num);
% y = zeros(M,2*ite_num);
% x_noDPD = zeros(M,2*ite_num);
% y_noDPD = zeros(M,2*ite_num);
% error_vec = zeros(M,ite_num);
% phi2 = zeros(M,4);
% phi2_noDPD = zeros(M,4);
% phi = zeros(2*(M+M^3),1);
    
% for nn = 2:ite_num
% %         s = s;
%     phi(1:M,1) = s(:,nn);
% %         for mm=1:M
% %             phi(mm*M+1:(mm+1)*M,1) = [s(mm)*abs(s).^2];
% % %             phi(:,mm) = [s;s(1)*abs(s).^2;s(2)*abs(s).^2];
% %         end
%     phi(M+1:M+M^3,1) = kron(kron(s(:,nn),s(:,nn)),conj(s(:,nn)));
%     phi(M+M^3+1:M+M^3+M,1) = s(:,nn-1);
%     phi(M+M^3+M+1:2*(M+M^3)) = kron(kron(s(:,nn-1),s(:,nn-1)),conj(s(:,nn-1)));
% 
%     for mm=1:M
%         x(mm,nn) = w(:,mm)'*phi;
%     end
%     if norm(x(:,nn))^2/M>2
%         w=zeros(2*(M+M^3),M);
%         for mm=1:M
%             w(mm,mm) = 1;
%         end
%     end
% 
% 
%     if EZmodel
%         % equivalent NL model in beamspace
%         phi2(:,1) = x(:,nn);
%         phi2(:,3) = x(:,nn-1);
%         for mm=1:M
%             phi2(mm,2) = x(mm,nn)*abs(x(mm,nn))^2+...
%                 2*x(mm,nn)*(sum(abs(x(:,nn)).^2)-abs(x(mm,nn))^2); 
% 
%             phi2(mm,4) = x(mm,nn-1)*abs(x(mm,nn-1))^2+...
%                 2*x(mm,nn-1)*(sum(abs(x(:,nn-1)).^2)-abs(x(mm,nn-1))^2); 
%         end
%         y(:,nn) = phi2*beta;%+noise(:,nn);
%     else
%         % Analog beamforming, NL PA, and combining
%         x_PA_ip = Frf*x;
%         x_PA_op = beta1*x_PA_ip+beta3*x_PA_ip.*abs(x_PA_ip).^2;
%         y(:,nn) = Frf'*x_PA_op/Nt;%+noise(:,nn);
%     end
% 
%     % calculate error vector and MSE
%     error_vec(:,nn) = s(:,nn)-y(:,nn);
%     MSE(nn,1) = norm(error_vec(:,nn))^2;
%     if MSE(nn,1)>100
%         stop = 1;
%     end
% 
%     % LMS updating based on partial derivative of r_k[n] to x_l[n-d]
%     % d = 0 part
%     for kk=1:M
%         for ll=1:M
%             if kk==ll
%                 h(ll,kk) = beta1_0+2*beta3_0*sum(abs(x(:,nn)).^2);
%             else
%                 h(ll,kk) = 2*beta3_0*x(ll,nn)*conj(x(kk,nn));
%             end
%         end
%     end
%     for mm=1:M
%         w(1:(M+M^3),mm) = w(1:(M+M^3),mm)+stepsize*(error_vec(:,nn))'*h(mm,:).'...
%             *(phi(1:(M+M^3)));
%     end
% 
%     % d = 1 part
%     for kk=1:M
%         for ll=1:M
%             if kk==ll
%                 h(ll,kk) = beta1_0+2*beta3_0*sum(abs(x(:,nn-1)).^2);
%             else
%                 h(ll,kk) = 2*beta3_0*x(ll,nn-1)*conj(x(kk,nn-1));
%             end
%         end
%     end
%     for mm=1:M
%         w((M+M^3+1):2*(M+M^3),mm) = w((M+M^3+1):2*(M+M^3),mm)+stepsize*(error_vec(:,nn))'*h(mm,:).'...
%             *(phi((M+M^3+1):2*(M+M^3)));
%     end
% end
% % figure;
% % plot(MSE)
%% DPD using weights
runtimes = 5e1;
ite_num = 1e3;
% C_withDPD = zeros(ite_num,runtimes);
C_noDPD = zeros(ite_num,runtimes);
C_ideal = zeros(ite_num,runtimes);

for runindex = 1:runtimes
    if mod(runindex,50)==0
        SIM_INFO = ['simulating',num2str(runindex/runtimes)];
        display(SIM_INFO);
    end
    
    
    P = ones(M,1)/M;
    noise_pow = ones(M,1)/SNR;
    %initialization of vector, matrix
    s = zeros(M,ite_num);
    x = zeros(M,1);

    
    for mm=1:M
        s(mm,:) = ones(1,ite_num)*sqrt(P(mm));
        s_raw = ((randn(1,ite_num))+1j*(randn(1,ite_num)))/sqrt(2);
%         s_raw = exp(1j*rand(1,ite_num)*2*pi)/sqrt(2);
%         s_raw = rand(1,ite_num)*6-3;
        PAPR = 2.5;
        index = find(abs(s_raw).^2<(PAPR));
        s(mm,index) = s_raw(index)*sqrt(P(mm));
    end

    x = zeros(M,ite_num);
    y = zeros(M,2*ite_num);
    x_noDPD = zeros(M,ite_num);
    y_noDPD = zeros(M,ite_num);
    error_vec = zeros(M,ite_num);
    phi2 = zeros(M,4);
    phi2_noDPD = zeros(M,4);
    phi = zeros(2*(M+M^3),1);
    
    
    for nn = 2:ite_num
%         phi(1:M,1) = s(:,nn);
%         phi(M+1:M+M^3,1) = kron(kron(s(:,nn),s(:,nn)),conj(s(:,nn)));
%         phi(M+M^3+1:M+M^3+M,1) = s(:,nn-1);
%         phi(M+M^3+M+1:2*(M+M^3)) = kron(kron(s(:,nn-1),s(:,nn-1)),conj(s(:,nn-1)));
%         % with DPD
%         for mm=1:M
%             x(mm,nn) = w(:,mm)'*phi;
%         end
%         phi2(:,1) = x(:,nn);
%         phi2(:,3) = x(:,nn-1);
%         for mm=1:M
%             phi2(mm,2) = x(mm,nn)*abs(x(mm,nn))^2+...
%                 2*x(mm,nn)*(sum(abs(x(:,nn)).^2)-abs(x(mm,nn))^2); 
% 
%             phi2(mm,4) = x(mm,nn-1)*abs(x(mm,nn-1))^2+...
%                 2*x(mm,nn-1)*(sum(abs(x(:,nn-1)).^2)-abs(x(mm,nn-1))^2); 
%         end
%         y(:,nn) = phi2*beta;
        
        % no DPD
        phi2_noDPD(:,1) = s(:,nn);
        phi2_noDPD(:,3) = s(:,nn-1);
        for mm=1:M
            phi2_noDPD(mm,2) = s(mm,nn)*abs(s(mm,nn))^2+...
                2*s(mm,nn)*(sum(abs(s(:,nn)).^2)-abs(s(mm,nn))^2); 

            phi2_noDPD(mm,4) = s(mm,nn-1)*abs(s(mm,nn-1))^2+...
                2*s(mm,nn-1)*(sum(abs(s(:,nn-1)).^2)-abs(s(mm,nn-1))^2); 
        end
        y_noDPD(:,nn) = phi2_noDPD*beta;
        for mm=1:M
%         int_pow_DPD = sum(abs(y(mm,nn)-s(mm,nn)).^2);
        int_pow_noDPD = sum(abs(y_noDPD(mm,nn)-s(mm,nn)).^2);

%         C_withDPD(nn,runindex) = C_withDPD(nn,runindex)+log2(1+P(mm)/(noise_pow(mm)+int_pow_DPD));
        C_noDPD(nn,runindex) = C_noDPD(nn,runindex)+log2(1+P(mm)/(noise_pow(mm)+int_pow_noDPD));
        C_ideal(nn,runindex) = C_ideal(nn,runindex)+log2(1+P(mm)/noise_pow(mm));
        end
        
    end
    
end
%%
% results(1,beam_index) = mean(mean(C_withDPD(10:end,:)));
results(2,beam_index) = mean(mean(C_noDPD(10:end,:)));
results(3,beam_index) = mean(mean(C_ideal(10:end,:)));
end
%%
% figure(1)
% plot(y,'o');hold on
% axis([-10,10,-10,10])
% grid on
%%
% figure(3)
% semilogy(mean(MSE,2));hold on
% grid on
% xlabel('Iteration Number')
% ylabel('MSE')
%%
% M = 2
% results(:,1)=[7.8053,5.3345,11.3449]';
% M = 3
% results(:,2)=[10.1750,7.2000,15.3046]';
% M = 4
% results(:,3)=[13.7431,9.0017,18.8018]';
% M = 5
% results(:,4)=[16.6748,10.7308,21.9616]';
% M = 6
% results(:,5)=[19.2479,12.4076,24.8577]';
figure(99)
plot(M_pool,results,'-o');hold on