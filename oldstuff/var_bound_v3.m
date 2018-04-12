%  05/22/2015
%  this script test variance of estimation by Monte Carlo Simulation
%  toghther with my theoretically upper bound
clear all;clc;

runtimes_sim = 2e4;
rrcPulseShape=1;
Target=1;

blk_num = 31;
M = 20;
N = 31;
SampleNumberAve = 30;
stat_num = 2e3;
P = M*N;
L = M*SampleNumberAve*stat_num;

%% sweeping one particular blk's bandwidth
R_range = fliplr([4,8,16,32,64,128,256,512,1024]);

for Rindex=1:length(R_range)
    Rindex
    SampleperSymbol=R_range(Rindex);
    SampleperFrame=M*SampleperSymbol;

    hdesign  = fdesign.pulseshaping(SampleperSymbol,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    rrc=hpulse.Numerator./sqrt(mean(hpulse.Numerator.^2));%./sqrt(length(hpulse.Numerator));
    for ii=1:length(rrc)
        h(ii)=rrc(ii:end)*rrc(1:end-ii+1)'/(length(rrc)-ii+1);
    end
    %  Define length of BLK we want.
    stat_num=40;
    L=SampleperFrame*stat_num;
    clearvars sig_bb BLK_r BLK_i BLK_sweep
    if rrcPulseShape
        for ii=1:50
            symbols=1500;
            hmod = modem.qammod('M', 4, 'InputType', 'integer');
            data = randi(4,symbols,1)-1;
            data = modulate(hmod, data);
            data = upsample(data,SampleperSymbol);
            temp_data = conv(data,hpulse.Numerator);
            sig_bb(:,ii) = temp_data(end-(symbols-5)*SampleperSymbol+1:end-5*SampleperSymbol+1);

            BLK_r(:,ii)=real(sig_bb(:,ii))./sqrt(mean(real(sig_bb(:,ii)).^2));
            BLK_i(:,ii)=imag(sig_bb(:,ii))./sqrt(mean(imag(sig_bb(:,ii)).^2));
            BLK_sweep(:,ii)=(BLK_r(:,ii)+BLK_i(:,ii))/sqrt(2);
        end
    else
        clearvars sig_bb_r sig_bb_i BLK_sweep
        for ii=1:50
            sig_bb_r(:,ii) = kron(randi(2,1500,1)*2-3,ones(SampleperSymbol,1));
            sig_bb_i(:,ii) = kron(randi(2,1500,1)*2-3,ones(SampleperSymbol,1));
            BLK_sweep(:,ii) = (sig_bb_r(:,ii)+sig_bb_i(:,ii))/sqrt(2);
        end    
    end

% MC estimation variance
    for runindex=1:runtimes_sim
        CAL = randi(2,1*N,blk_num)*2-3;

        start_point = randi(1400*SampleperSymbol,blk_num);
        sig_cache = BLK_sweep(start_point(ii):start_point(ii)+N*1-1,randi(50));
        
        for ii=1:blk_num
            if Target
                result_sim(runindex) = tagging_v3(sig_cache,ones(N*1,1),N,1);
            else
                result_sim(runindex) = tagging_v3(sig_cache.*CAL(:,1),CAL(:,2),N,1);
            end
        end
    end
    if rrcPulseShape
        mus2_anal(Rindex)=anal_mean_fun_v2(SampleperSymbol,N,h,'t');
    else
        mus2_anal(Rindex)=anal_mean_fun(SampleperSymbol,N,'t');
    end
    var_sim(Rindex) = var(result_sim);
    mean_sim(Rindex) = mean(result_sim);
end
%% plot
if Target
    ydata = var_sim-2*mean_sim.^2;
    xdata=4./R_range;
    figure
    plot_setting
    plot(xdata(1:1:end),ydata(1:1:end),'bo-');hold on
    plot(4./R_range,-2*mus2_anal.^2,'b--');hold on
    plot(4./R_range,-3*mus2_anal.^2+2.233,'b--');hold on
    grid on
    xlabel('Interferer Bandwidth');
    ylabel('l^th Diagonal Component of  \Sigma');
    legend('Sim.(RRC)','Lower Approx.(RRC)','Upper Approx.(RRC)')
else
    ydata = var_sim-2*mean_sim.^2;
    xdata=4./R_range;
    figure
    plot_setting
    plot(xdata(1:1:end),ydata(1:1:end),'bo-');hold on
    plot(4./R_range,ones(1,length(R_range))*3*(0.972-1)/N^2,'b--');hold on
    plot(4./R_range,ones(1,length(R_range))*3*(1.986-1)/N^2,'b--');hold on
    grid on
    xlabel('Interferer Bandwidth');
    ylabel('i^th Diagonal Component of  \Sigma');
    legend('Sim.(RRC)','Lower Approx.(RRC)','Upper Approx.(RRC)')
end
