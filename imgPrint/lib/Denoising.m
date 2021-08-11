function [DenoisingSignal] = Denoising(ecgData)
    DenoisingSignal = [];
    for i = 1:size(ecgData,2)
        level=8; wavename='coif5';
        [C,L]=wavedec(ecgData(:,i),level,wavename);
        A8=appcoef(C,L,wavename,8);
        D1=detcoef(C,L,1);
        D2=detcoef(C,L,2);
        D3=detcoef(C,L,3);
        D4=detcoef(C,L,4);
        D5=detcoef(C,L,5);
        D6=detcoef(C,L,6);
        D7=detcoef(C,L,7);
        D8=detcoef(C,L,8);
        A8=zeros(length(A8),1); %ȥ������Ư��,8�ӵ��l��Ϣ
        RA7=idwt(A8,D8,wavename);
        RA6=idwt(RA7(1:length(D7)),D7,wavename);
        RA5=idwt(RA6(1:length(D6)),D6,wavename);
        RA4=idwt(RA5(1:length(D5)),D5,wavename);
        RA3=idwt(RA4(1:length(D4)),D4,wavename);
        RA2=idwt(RA3(1:length(D3)),D3,wavename);
        D2=zeros(length(D2),1); %ȥ�����l����2�Ӹ��l��
        RA1=idwt(RA2(1:length(D2)),D2,wavename);
        D1=zeros(length(D1),1);%ȥ�����l����1�Ӹ��l��
        DenoisingSignal(:,i)=idwt(RA1(1:length(D1)),D1,wavename); %��ȷ���᲻����ʧ��Ϣ
    end
end