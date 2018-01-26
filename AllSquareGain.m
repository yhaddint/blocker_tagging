% 03/04/2015
% Check Suppression Gain when all waveform is square. MC v.s. Theoretical
% Done! Now I will compare with rrc waveform

clear;
clc;
clf
warning off


% Independent repeating experiment.
runtimes=1e4;

N_range=[15 31 63];
N_num = length(N_range);
M=50;

    % BLK
R_range = 10:4:120;
for Rindex = 1:length(R_range)
        
    SampleperSymbol=R_range(Rindex);
    SampleperFrame=M*SampleperSymbol;
    
    % theoretical value
    % First to formulate B matrix
%     B = zeros(SampleperSymbol,SampleperSymbol);
%     for rr=1:SampleperSymbol
%         value = (SampleperSymbol-rr+1)/(SampleperSymbol);
%         B = B + diag(ones(1,SampleperSymbol-rr+1).*value,rr-1);
%         if rr > 1
%             B = B + diag(ones(1,SampleperSymbol-rr+1).*value,-(rr-1));
%         end
%     end
%     
    %  How many BLK we want to tag simutaneously
    tag_num=20;

    %  How many active BLK existing. To get reasonable number of different
    %  envelope, we generate 2x of them.
    BLK_num=10;
    BLKcan_num=20;

    %  Define length of BLK we want.
    stat_num=40;
    L=SampleperFrame*stat_num;
    clearvars sig_bb_r sig_bb_i
    for ii=1:BLKcan_num
        sig_bb_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleperSymbol,1));
        sig_bb_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleperSymbol,1));
    end

    for nindex=1:N_num

        %  parameter M and N setting. N is also length of CAL Sequences
        N=N_range(nindex);
        
        % BLK correlation matrix
        B_corr = zeros(N,N);
        for rr=1:N
            if rr <= SampleperSymbol
                value = (SampleperSymbol-rr+1)/(SampleperSymbol);
                B_corr = B_corr+diag(ones(1,N-rr+1).*value,rr-1);
                if rr > 1
                    B_corr = B_corr + diag(ones(1,N-rr+1).*value,-(rr-1));
                end
            end
        end
        
        % CAL correlation matrix
        C_corr = ones(N,N).*(-1/N);
        C_corr = C_corr + diag(ones(1,N)*(N+1)/N,0);
        
        Gain_theory(nindex,Rindex) = 10*log10(sum(sum(B_corr))/trace(B_corr*C_corr));


        % calibration sequences
        cal0=PN_generator_downsample(N,N,1);
        CAL=LowerRate(cal0,1,N);

        %% Interference Generation
        clearvars BLK
        for runindex=1:runtimes
            %  permutation of calibration sequences
            CAL_index=randi([2,N]);

            tempindex=randi(BLKcan_num);

            Starindex=randi(L-SampleperFrame-1);
            
            BLK(:,1)=sig_bb_r(Starindex:Starindex+N-1,tempindex)+sig_bb_i(Starindex:Starindex+N-1,tempindex);

            Power_Target(runindex,nindex,Rindex) = correlation_and_pow(BLK(:,1),ones(size(BLK(:,1))),N,1);
            Power_NonTarget(runindex,nindex,Rindex) = correlation_and_pow(BLK(:,1).*CAL(:,1),CAL(:,CAL_index),N,1);
        end
    Gain(nindex,Rindex) = 10*log10(sum(squeeze(Power_Target(:,nindex,Rindex)))/sum(squeeze(Power_NonTarget(:,nindex,Rindex))));
    end
end
%% plot
figure;plot_setting();
plot(R_range,Gain_theory(1,:),'b');hold on
plot(R_range,Gain_theory(2,:),'r');hold on
plot(R_range,Gain_theory(3,:),'g');hold on
plot(R_range,Gain(1,:),'bx');hold on
plot(R_range,Gain(2,:),'rx');hold on
plot(R_range,Gain(3,:),'gx');hold on
legend('N = 15 (Anal.)','N = 31 (Anal.)','N = 63 (Anal.)','N = 15 (Sim.)','N = 31 (Sim.)','N = 63 (Sim.)');
grid on
xlabel('BW ratio');
ylabel('Suppression (dB)')
