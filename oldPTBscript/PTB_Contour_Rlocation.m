clc; clear;close all;
%%
global sfreq points
fprintf(1,'PDB DATA LOADING');
%------ SPECIFY DATA ------------------------------------------------------
stringname = "s0260lre"; %�ļ���
patientNumber = "080";
points= 10000;%����Ѷ�ŵ������
PATH= pwd+"\ptbdb\patient"+patientNumber+"\"; % ����·��
HEADERFILE= strcat(stringname,'.hea');     % header-file in text format
ATRFILE= strcat(stringname,'.atr');        % attributes-file in binary format
DATAFILE=strcat(stringname,'.dat');        % data-file
XYZFILE=strcat(stringname,'.xyz');        % data-file
leadname=["I","II","III","avr","avl","avf","v1","v2","v3","v4","v5","v6","vx","vy","yz"];
%------ LOAD HEADER DATA --------------------------------------------------
fprintf(1,'\nWORKING ON %s ...\n', HEADERFILE);
signalh= fullfile(PATH, HEADERFILE);
fid1=fopen(signalh,'r');
z= fgetl(fid1);
A= sscanf(z, '%*s %d %d %d',[1,3]);
sfreq=A(2);   % sample rate of dataȡ��Ƶ��
clear A;
for k=1:15
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
    dformat(k)= A(1);           %#ok<*SAGROW> % format; here only 212 is allowed
    gain(k)= A(2);              % number of integers per mV
    bitres(k)= A(3);            % bitresolution
    zerovalue(k)= A(4);         % integer value of ECG zero point
    firstvalue(k)= A(5);        % first integer value of signal (to test for errors)
end
fclose(fid1);
clear A;
M = loadfile(fullfile(PATH, DATAFILE),12);
XYZ = loadfile(fullfile(PATH, XYZFILE),3);
M=[M XYZ];
for i = 1:15
    M(:,i)=(M(:,i)-zerovalue(i))/gain(i);
end
TIME = (0:(points-1))/sfreq;
%-------------------------------------------------
%{
figure(1)
[~,c]=contourf(1:points,1:15,M(1:points,:)',100);
set(c,'LineColor','none')
title("patient"+patientNumber+" "+stringname)
colormap(jet)
%}
for lead= 1:15
    %clearvars -except M points sfreq stringname leadname lead
    ecgdata=M(:,lead);
    for i = 2:length(ecgdata)-1   %��ֵ�˲�
        ecgdata(i)=(ecgdata(i-1)+ecgdata(i)+ecgdata(i+1))/3;
    end
    %{
    %���fig����
    if lead>1
        all_figs = findobj(0, 'type', 'figure');
        want_figs = 16;
        for fignu = 1:fix((lead-2)/3)
            want_figs = [want_figs 16+fignu];
        end
        clean_figs = setdiff(all_figs,want_figs);
        delete(setdiff(all_figs,want_figs));
    end
    %}
    
    %%
    %%%%%%%%%%%%%%%%%%%ȥ����Ѷ�����ͻ���Ư��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    level=8; wavename='coif5';
    % level=8; wavename='db5';
    %figure(2);
    %plot(ecgdata(1:points));grid on ;axis tight;axis([1,points,-2,5]);
    %title('ԭʼECGѶ��');
    %%%%%%%%%%����С��ת��8��
    [C,L]=wavedec(ecgdata,level,wavename);
    %%%%%%%��ȡС���߶�ϵ����
    A1=appcoef(C,L,wavename,1);
    A2=appcoef(C,L,wavename,2);
    A3=appcoef(C,L,wavename,3);
    A4=appcoef(C,L,wavename,4);
    A5=appcoef(C,L,wavename,5);
    A6=appcoef(C,L,wavename,6);
    A7=appcoef(C,L,wavename,7);
    A8=appcoef(C,L,wavename,8);
    %%%%%%%��ȡϸ��ϵ��
    D1=detcoef(C,L,1);
    D2=detcoef(C,L,2);
    D3=detcoef(C,L,3);
    D4=detcoef(C,L,4);
    D5=detcoef(C,L,5);
    D6=detcoef(C,L,6);
    D7=detcoef(C,L,7);
    D8=detcoef(C,L,8);
    %%%%%%%%%%%%�ع�
    A8=zeros(length(A8),1); %ȥ������Ư��,8���Ƶ��Ϣ
    RA7=idwt(A8,D8,wavename);
    RA6=idwt(RA7(1:length(D7)),D7,wavename);
    RA5=idwt(RA6(1:length(D6)),D6,wavename);
    RA4=idwt(RA5(1:length(D5)),D5,wavename);
    RA3=idwt(RA4(1:length(D4)),D4,wavename);
    RA2=idwt(RA3(1:length(D3)),D3,wavename);
    D2=zeros(length(D2),1); %ȥ����Ƶ������2���Ƶ����
    RA1=idwt(RA2(1:length(D2)),D2,wavename);
    D1=zeros(length(D1),1);%ȥ����Ƶ������1���Ƶ����
    DenoisingSignal=idwt(RA1,D1,wavename);
    %figure(3);
    %plot(DenoisingSignal);
    %title('ȥ����ƵѶ�������Ʈ�Ƶ�ECGѶ��'); grid on; axis tight;axis([1,points,-2,5]);
    M(:,lead) =  DenoisingSignal;
    clear ecgdata;
end

%clearvars -except M points sfreq stringname leadname lead patientNumber
figure(777)
[M1,c]=contour(1:points,1:15,M(1:points,:)',100);
set(c,'LineColor','none')
set(c,'Fill','on')
set(c,'Levellist',-2:0.01:2)
title("patient"+patientNumber+" "+stringname)
%colormap(jet)
Rloc = [0 0 0 0];    
FindingR
Rloc = mean(rmoutliers(Rloc),1);%%ɾ����Ⱥֵ����ƽ��
Rloc (:,all(Rloc == 0, 1)) = [];%%ɾ����ֵ
line(repelem(Rloc,15,1),1:15,"Color","white","LineWidth",1,'LineStyle','- -');
%%
function data=loadfile(filepath,nosig)
global points
fid2=fopen(filepath,'r');
data= fread(fid2, [nosig, points], 'uint16')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);                                 % ������3�У�ÿ��8λ��= 2 * 12λ
for i=1:nosig
    PR = bitand(data(:,i),32768);
    data(:,i)=bitand(data(:,i),32767)-PR;
end
end