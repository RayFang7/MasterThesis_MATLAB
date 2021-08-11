function [DenoisingSignal] = Denoising2(ecgData)
%     for i = 1:size(ecgData,1)-6
%         ecgData(i+4,:) = sum(ecgData(i:i+6,:))/7;
%     end
    DenoisingSignal = [];
    for i = 1:size(ecgData,2)
        level=9; wavename='db6';
        [C,L]=wavedec(ecgData(:,i),level,wavename);
        A9=appcoef(C,L,wavename,9);
        D1=detcoef(C,L,1);
        D2=detcoef(C,L,2);
        D3=detcoef(C,L,3);
        D4=detcoef(C,L,4);
        D5=detcoef(C,L,5);
        D6=detcoef(C,L,6);
        D7=detcoef(C,L,7);
        D8=detcoef(C,L,8);
        D9=detcoef(C,L,9);
        A9=zeros(length(A9),1);
        RA8=idwt(A9,D9,wavename);
        RA7=idwt(RA8(1:length(D8)),D8,wavename);
        RA6=idwt(RA7(1:length(D7)),D7,wavename);
        RA5=idwt(RA6(1:length(D6)),D6,wavename);
        RA4=idwt(RA5(1:length(D5)),D5,wavename);
        RA3=idwt(RA4(1:length(D4)),D4,wavename);
        D3=zeros(length(D3),1); %去除高頻噪聲，2層高頻噪聲
        RA2=idwt(RA3(1:length(D3)),D3,wavename);
        D2=zeros(length(D2),1); %去除高頻噪聲，2層高頻噪聲
        RA1=idwt(RA2(1:length(D2)),D2,wavename);
        D1=zeros(length(D1),1);%去除高頻噪聲，1層高頻噪聲
        DenoisingSignal(:,i)=idwt(RA1(1:length(D1)),D1,wavename); %不确定会不会损失信息
    end
end