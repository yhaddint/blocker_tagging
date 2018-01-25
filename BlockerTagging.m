% 03/25/2014
%
% This code is to generate QSKP signal by comm_toolbox and verify 
% certian CLT approximation valid or not

% 03/27/2014
% PN sequences with characterize polynomial have been put into use

% 04/11/2014
% All noncancelling and cancelling terms are considered
% Found there is '.^2' in sig_use, which was wrong. Fixed
% Surprisingly enough, statistics of noncancelling terms is
% gaussian now

clear;clc;clf;close all
warning off
symbols=3e4;
upsam=8;
Npool=[20 40 60 80];
Mpool=[1 2 4 8 12 16];
% Npool=[100];
% Mpool=[4];


needmoredata=1;
height=2e3;
rslt_matrix=zeros(length(Npool),length(Mpool),height);
mark=zeros(length(Npool),length(Mpool));

for runtime=1:5  
    for i=1:1
    clear data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    data = conv(data,hpulse.Numerator);

    sig_bb(:,i)=data;
    sig_bb(:,i)=sig_bb(:,i)/(sqrt(mean(abs(sig_bb(:,i)).^2)));
    end





    %% distribution of signal_r^2*cal
    L=2e5;
    sig_use=real(sig_bb(1:L,1));
    sig_power=mean(real(sig_bb(:,1)).^2);
    cal=PNgenerator_v1(L);
    sig_pn=sig_use.*cal;
    
    
    for nn=1:length(Npool)
    for mm=1:length(Mpool)
        clc
        display(['Run number',num2str(runtime)]);
        display(['N = ',num2str(Npool(nn)) ', M = ',num2str(Mpool(mm))]);
    clear rslt_cancel rslt_noncancel rslt_PN ii i
    
    N=Npool(nn);
    M=Mpool(mm);

    for ii=1:fix(L/N/2)
        rslt_cancel(ii)=mean(sig_pn((ii-1)*N+1:ii*N));
        rslt_noncancel(ii)=mean(sig_use((ii-1)*N+1:ii*N));
        rslt_PN(ii)=mean(cal((ii-1)*N+1:ii*N));
    end
    %rslt3=abs(rslt_cancel+rslt_noncancel+rslt_PN);
    rslt3=rslt_noncancel;
    
    l_need=fix(length(rslt3)/M);
    if M==1
        temp=rslt3;
    else
        temp=mean(reshape(rslt3(1:M*l_need),M,l_need));
    end
    
    
    if mark(nn,mm)<height
    % it means we should update new data into it
        if mark(nn,mm)+length(temp)<=height
            rslt_matrix(nn,mm,mark(nn,mm)+1:mark(nn,mm)+length(temp))=temp;
            mark(nn,mm)=mark(nn,mm)+length(temp);
        else
            rslt_matrix(nn,mm,mark(nn,mm)+1:height)=temp(1:height-mark(nn,mm));
            mark(nn,mm)=height;
        end
    end
    %%

    % figure
    % hist(rslt2,20)
    % title('sequencially pick')


    %%
    % clear ii
    % for ii=1:fix(N/M)
    %     rslt4(ii)=mean(cal((ii-1)*M+1:ii*M).^3);
    % end

    % %% sequencially pick illustration
    % figure
    % plot(real(sig_bb(1:500,1)),'c','linewidth',2);
    % hold
    % plot(100:150,real(sig_bb(100:150,1)),'kx');
    % grid
    % 
    % 
    % %% randomly pick illustration
    % index=randi(500,1,50);
    % figure
    % plot(real(sig_bb(1:500,1)),'c','linewidth',2);
    % hold
    % plot(sort(index),real(sig_bb(sort(index))),'kx');
    % grid


    %[H,gaussfit(mm,nn,:)]=chi2gof(rslt_matrix(mm,nn,:));
    end
    end
end
%% surf of GOF
for nn=1:length(Npool)
    for mm=1:length(Mpool)
        temp=squeeze(rslt_matrix(nn,mm,:));
        G_mean(nn,mm)=mean(temp(find(temp)));
        G_var(nn,mm)=mean(temp(find(temp)).^2);
        [H rslt_hist(nn,mm)]=chi2gof(temp(find(temp)));
    end
end
rslt_hist
surf(Mpool,Npool,rslt_hist);

%% hist of noncancelling term
figure
for ii=1:4
subplot(2,2,ii)
hist(squeeze(rslt_matrix(ii,1,:)))
title(['N=',num2str(Npool(ii))]);
end