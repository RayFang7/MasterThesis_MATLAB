clc; clear;close all;
%%
global sfreq points
fprintf(1,'PDB DATA LOADING');
%------ SPECIFY DATA ------------------------------------------------------
stringname = "s0260lre"; %�ļ���
patientNumber = "080";
points= 38400;%����Ѷ�ŵ������
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
for te= 1:15
    clearvars -except M te points sfreq stringname leadname
    ecgdata=M(:,te);
    if te>1
        all_figs = findobj(0, 'type', 'figure');
        want_figs = 16;
        for fignu = 1:fix((te-2)/3)
            want_figs = [want_figs 16+fignu];
        end
        clean_figs = setdiff(all_figs,want_figs);
        delete(setdiff(all_figs,want_figs));
    end
    %%
    %%%%%%%%%%%%%%%%%%%ȥ����Ѷ�����ͻ���Ư��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    level=8; wavename='coif5';
    % level=8; wavename='db5';
    
    figure(2);
    plot(ecgdata(1:points));grid on ;axis tight;axis([1,points,-2,5]);
    title('ԭʼECGѶ��');
    for i = 2:length(ecgdata)-1   %��ֵ�˲�
        ecgdata(i)=(ecgdata(i-1)+ecgdata(i)+ecgdata(i+1))/3;
    end
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
    
    D1=detcoef(C,L,1);
    D2=detcoef(C,L,2);
    D3=detcoef(C,L,3);
    D4=detcoef(C,L,4);
    D5=detcoef(C,L,5);
    D6=detcoef(C,L,6);
    D7=detcoef(C,L,7);
    D8=detcoef(C,L,8);
    
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
    plot(DenoisingSignal);
    title('ȥ����ƵѶ�������Ʈ�Ƶ�ECGѶ��'); grid on; axis tight;axis([1,points,-2,5]);
    clear ecgdata;
    %%%%%%%%%%%%%
    % figure(5000)
    % plot(A8)
    % hold on
    % plot(initialA8)
    % hold off
    %%%%%%%%%%%%%%
    % <pre name="code" class="cpp">
    level=4;    sr=sfreq;
    %����ECG�ź�
    %load ecgdata.mat;
    %load ECGsignalM1.mat;
    %load Rsignal.mat
    mydata = DenoisingSignal;
    ecgdata=mydata';
    swa=zeros(4,points);%�����ò��Ϣ    ����ͬ����С�������
    swd=zeros(4,points);%����ϸ����Ϣ    ����ͬ����С�������
    signal=ecgdata(0*points+1:1*points); %ȡ���ź�
    
    %��С��ϵ���ͳ߶�ϵ��
    %��ͨ�˲��� 1/4 3/4 3/4 1/4
    %��ͨ�˲��� -1/4 -3/4 3/4 1/4
    %��������С��
    %%%%%%%%%%%%%%%%%%%%% ��ɢС������---����С�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:points-3
        swa(1,i+3)=1/4*signal(i+3-2^0*0)+3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
        swd(1,i+3)=-1/4*signal(i+3-2^0*0)-3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
    end
    j=2;
    while j<=level
        for i=1:points-24
            swa(j,i+24)=1/4*swa(j-1,i+24-2^(j-1)*0)+3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
            swd(j,i+24)=-1/4*swa(j-1,i+24-2^(j-1)*0)-3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
        end
        j=j+1;
    end
    %%%%%%%%%%%%%%%%%%%%%% ��ɢС������---����С�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ECGѶ�ż�����j=1,2,3,4�߶��µĳ߶�ϵ����С��ϵ��
    %{
figure(11);
subplot(level,1,1); plot(ecgdata(1:points)); grid on;axis tight;
title('ECGѶ�ż�����j=1,2,3,4�߶��µĳ߶�ϵ����С��ϵ��');
    %}
    for i=1:level
        subplot(level+1,2,2*(i)+1);
        plot(swa(i,:)); axis tight;grid on;xlabel('time');
        ylabel(strcat('a   ',num2str(i)));
        subplot(level+1,2,2*(i)+2);
        plot(swd(i,:)); axis tight;grid on;
        ylabel(strcat('d   ',num2str(i)));
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%% P��Q��R��S��T���ͱ߽��ж� %%%%%%%%%%%%%%%%%%%%%%%%%%
    %************************** ����������ֵ **********************************%
    ddw=zeros(size(swd));    %����ͬ����С�������
    pddw=ddw;
    nddw=ddw;
    %%%%%%%%%%%%%%%%%%%%%%%%% ����ֵ�ļ�ֵ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %С��ϵ���Ĵ���0��?
    posw=swd.*(swd>0);
    %б�ʴ���0
    pdw=((posw(:,1:points-1)-posw(:,2:points))<0);
    %������ֵ��
    pddw(:,2:points-1)=((pdw(:,1:points-2)-pdw(:,2:points-1))>0);
    %%%%%%%%%%%%%%%%%%%%%%%%%% �Ҹ�ֵ�ļ�ֵ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %С��ϵ��С��0�ĵ�
    negw=swd.*(swd<0);
    ndw=((negw(:,1:points-1)-negw(:,2:points))>0);
    %������ֵ��
    nddw(:,2:points-1)=((ndw(:,1:points-2)-ndw(:,2:points-1))>0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %������
    ddw=pddw|nddw;
    ddw(:,1)=1;
    ddw(:,points)=1;
    %�����ֵ���ֵ,����������Ϊ0
    wpeak=ddw.*swd;
    wpeak(:,1)=wpeak(:,1)+1e-10;
    wpeak(:,points)=wpeak(:,points)+1e-10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % �������߶��¼�ֵ��
    %{
figure(13);
for i=1:level
   subplot(level,1,i);
   plot(wpeak(i,:)); axis tight;grid on;
ylabel(strcat('j=   ',num2str(i)));
end
subplot(4,1,1);
title('ECG�ź���j=1,2,3,4�߶��µ�С��ϵ����ģ����ֵ��');
    %}
    interva2=zeros(1,points);
    intervaqs=zeros(1,points);
    Mj1=wpeak(1,:);
    Mj3=wpeak(3,:);
    Mj4=wpeak(4,:);
    %�����߶�3��ֵ��
    figure(14);
    plot (Mj4);axis tight;
    %title('�߶�4��С��ϵ����ģ����ֵ��');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% ����Mj3���Mj4 %%%%%%%%
    posi=Mj4.*(Mj4>0);
    %��������ֵ��ƽ��
    thposi=(max(posi(1:round(points/4)))+max(posi(round(points/4):2*round(points/4)))+max(posi(2*round(points/4):3*round(points/4)))+max(posi(3*round(points/4):4*round(points/4))))/4;
    posi=(posi>thposi/3);
    nega=Mj4.*(Mj4<0);
    %�󸺼���ֵ��ƽ��
    thnega=(min(nega(1:round(points/4)))+min(nega(round(points/4):2*round(points/4)))+min(nega(2*round(points/4):3*round(points/4)))+min(nega(3*round(points/4):4*round(points/4))))/4;
    nega=-1*(nega<thnega/4);
    %�ҳ���0��
    interva=posi+nega;
    loca=find(interva);
    for i=1:length(loca)-1
        if abs(loca(i)-loca(i+1))<80
            diff(i)=interva(loca(i))-interva(loca(i+1));
        else
            diff(i)=0;
        end
    end
    %�ҳ���ֵ-������
    loca2=find(diff==-2);
    %������ֵ��
    interva2(loca(loca2(1:length(loca2))))=interva(loca(loca2(1:length(loca2))));
    %������ֵ��
    interva2(loca(loca2(1:length(loca2))+1))=interva(loca(loca2(1:length(loca2))+1));
    intervaqs(1:points-10)=interva2(11:points);
    countR=zeros(1,1);
    countQ=zeros(1,1);
    countS=zeros(1,1);
    mark1=0;
    mark2=0;
    mark3=0;
    i=1;
    j=1;
    Rnum=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% PQSRT���ͼ�ⷶΧ�趨 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %*********** ��������ֵ�Թ��㡣��R����ֵ��������QRS����㼰�յ� ***********%
    % j=1;
    % Rnumx=[];
    % Rnumy=[];
    while i<points
        if interva2(i)==-1
            mark1=i;
            i=i+1;
            while(i<points&interva2(i)==0)
                i=i+1;
            end
            mark2=i;
            %�󼫴�ֵ�ԵĹ����
            %%%% ����Mj3���Mj4 %%%%%%%%
            mark3= round((abs(Mj4(mark2))*mark1+mark2*abs(Mj4(mark1)))/(abs(Mj4(mark2))+abs(Mj4(mark1))));
            [Rtestx,Rtesty]=max(ecgdata(mark3-40:mark3));
            Rtestyy=mark3-40+Rtesty;
            movepoint=Rtesty;
            Qmove=mark1-movepoint;
            Smove=mark2-movepoint;
            %        [Rnumx(j),Rnumy(j)]=max(ecgdata(mark3-40:mark3));
            %        j=j+1;
            %R������ֵ��
            %        R_result(j)=mark3-15;           %Ϊ�� -10�� ���Գ�����
            R_result(j)=Rtestyy;
            %        countR(mark3-15)=1;              %-10�ĳ�-15���Ŵ�Χ����
            countR(Rtestyy)=1;
            %���QRS�����
            %        kqs=mark3-15;                    %-10�ĳ�-15���Ŵ�Χ����
            %        kqs=Rtestyy-5;
            kqs=Qmove;  %%%% test
            markq=0;
            while (kqs>1)&&( markq< 3)
                if Mj1(kqs)~=0
                    markq=markq+1;
                end
                kqs= kqs -1;
            end
            countQ(kqs)=-1;
            %���QRS���յ�
            %   kqs=mark3+5;                     %-10�ĳ�-5���Ŵ�Χ����
            %   kqs=Rtestyy+5;
            kqs=Smove;   %%% test
            marks=0;                         %�ٸĳ�+5����Ϊ���ⲡ��QRSʱ�䳤�Ƚϳ���
            while (kqs<points)&&( marks<3)   %��S����ȡ���ֵ������յ㷶Χ�������ڱȽ�Զ����
                if Mj1(kqs)~=0
                    marks=marks+1;
                end
                kqs= kqs+1;
            end
            countS(kqs)=-1;
            i=i+60;
            j=j+1;
            Rnum=Rnum+1;
        end
        i=i+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %************************ ɾ�����㣬����©��� **************************%
    num2=1;
    while(num2~=0)
        num2=0;
        %j=3,�����
        R=find(countR);
        %�������
        R_R=R(2:length(R))-R(1:length(R)-1);
        RRmean=mean(R_R);
        %������R�����С��0.4RRmeanʱ,ȥ��ֵС��R��
        for i=2:length(R)
            if (R(i)-R(i-1))<=0.4*RRmean
                num2=num2+1;
                if signal(R(i))>signal(R(i-1))
                    countR(R(i-1))=0;
                else
                    countR(R(i))=0;
                end
            end
        end
    end
    
    num1=2;
    while(num1>0)
        num1=num1-1;
        R=find(countR);
        R_R=R(2:length(R))-R(1:length(R)-1);
        RRmean=mean(R_R);
        %������R���������1.6RRmeanʱ,��С��ֵ,����һ�μ��R��
        for i=2:length(R)
            if (R(i)-R(i-1))>1.6*RRmean
                Mjadjust=wpeak(4,R(i-1)+80:R(i)-80);
                points2=(R(i)-80)-(R(i-1)+80)+1;
                %��������ֵ��
                adjustposi=Mjadjust.*(Mjadjust>0);
                adjustposi=(adjustposi>thposi/4);
                %�󸺼���ֵ��
                adjustnega=Mjadjust.*(Mjadjust<0);
                adjustnega=-1*(adjustnega<thnega/5);
                %������
                interva4=adjustposi+adjustnega;
                %�ҳ���0��
                loca3=find(interva4);
                diff2=interva4(loca3(1:length(loca3)-1))-interva4(loca3(2:length(loca3)));
                %�����м���ֵ��,�ҳ�����ֵ��
                loca4=find(diff2==-2);
                interva3=zeros(points2,1)';
                for j=1:length(loca4)
                    interva3(loca3(loca4(j)))=interva4(loca3(loca4(j)));
                    interva3(loca3(loca4(j)+1))=interva4(loca3(loca4(j)+1));
                end
                mark4=0;
                mark5=0;
                mark6=0;
                while j<points2
                    if interva3(j)==-1;
                        mark4=j;
                        j=j+1;
                        while(j<points2&interva3(j)==0)
                            j=j+1;
                        end
                        mark5=j;
                        %������
                        mark6= round((abs(Mjadjust(mark5))*mark4+mark5*abs(Mjadjust(mark4)))/(abs(Mjadjust(mark5))+abs(Mjadjust(mark4))));
                        countR(R(i-1)+80+mark6-10)=1;
                        j=j+60;
                    end
                    j=j+1;
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%����ԭͼ����ʾ�������
    %%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ��PT����
    %��R����ǰ�Ĳ����üӴ����������СΪ100��Ȼ����㴰���ڼ���С�ľ���
    % figure(20);
    % plot(Mj4);
    % title('j=4 ϸ��ϵ��');  hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%����ֱ����j=4ʱ��R������
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mj4posi=Mj4.*(Mj4>0);
    %��������ֵ��ƽ��
    Mj4thposi=(max(Mj4posi(1:round(points/4)))+max(Mj4posi(round(points/4):2*round(points/4)))+max(Mj4posi(2*round(points/4):3*round(points/4)))+max(Mj4posi(3*round(points/4):4*round(points/4))))/4;
    Mj4posi=(Mj4posi>Mj4thposi/3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mj4nega=Mj4.*(Mj4<0);
    %�󸺼���ֵ��ƽ��
    Mj4thnega=(min(Mj4nega(1:round(points/4)))+min(Mj4nega(round(points/4):2*round(points/4)))+min(Mj4nega(2*round(points/4):3*round(points/4)))+min(Mj4nega(3*round(points/4):4*round(points/4))))/4;
    Mj4nega=-1*(Mj4nega<Mj4thnega/4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Mj4interval=Mj4posi+Mj4nega;   %�����뼫С��ϳ�һ��
    Mj4local=find(Mj4interval);    %�ҵ������λ�ã����м�������ֵ����ԭʼѶ�ŵ��ĸ���
    Mj4interva2=zeros(1,points);   %������Ѷ�ŵ���һ���������
    for i=1:length(Mj4local)-1
        if abs(Mj4local(i)-Mj4local(i+1))<80      %ͬһ��ECG������
            Mj4diff(i)=Mj4interval(Mj4local(i))-Mj4interval(Mj4local(i+1));
        else
            Mj4diff(i)=0;
        end
    end
    %�ҳ���ֵ��
    Mj4local2=find(Mj4diff==-2);
    %������ֵ��
    Mj4interva2(Mj4local(Mj4local2(1:length(Mj4local2))))=Mj4interval(Mj4local(Mj4local2(1:length(Mj4local2))));
    %������ֵ��
    Mj4interva2(Mj4local(Mj4local2(1:length(Mj4local2))+1))=Mj4interval(Mj4local(Mj4local2(1:length(Mj4local2))+1));
    mark1=0;
    mark2=0;
    mark3=0;
    Mj4countR=zeros(1,1);
    Mj4countQ=zeros(1,1);
    Mj4countS=zeros(1,1);
    flag=0;
    while i<points
        if Mj4interva2(i)==-1
            mark1=i;
            i=i+1;
            while(i<points&Mj4interva2(i)==0)
                i=i+1;
            end
            mark2=i;
            %%�󼫴�ֵ�ԵĹ����,��R4�м�ֵ֮���������R�㡣
            
            mark3= round((abs(Mj4(mark2))*mark1+mark2*abs(Mj4(mark1)))/(abs(Mj4(mark2))+abs(Mj4(mark1))));
            Mj4countR(mark3)=1;
            Mj4countQ(mark1)=-1;
            Mj4countS(mark2)=-1;
            flag=1;
        end
        if flag==1
            i=i+200;
            flag=0;
        else
            i=i+1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%�ҵ�MJ4��QRS�������ȱ�ٶ�R���©����Ͷ��ڼ�⡣�Ȳ�ȥϸ���ˡ�
    %%%%%
    %%%%%�Գ߶�4��R���ⲻ���ã���Ҫ�Ľ��ĵط�
    %%%%%%
    figure(200);
    plot(Mj4);
    title('�߶�4��R��С��ϵ����ǵ�');
    hold on;
    plot(ecgdata,'m');
    plot(Mj4countR,'r');
    plot(Mj4countQ,'g');
    plot(Mj4countS,'g');axis tight;
    % plot(ecgdata,'m');
    % plot(Mj4interva2);
    % plot(countR,'r');
    % plot(countQ,'g');
    % plot(countS,'g');
    axis tight;
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Mj4������ҵ� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Mj4countR.Q.S�ĳ�countR.Q.S%%%%%%%%%%%%%
    Rlocated=find(Mj4countR);
    Qlocated=find(Mj4countQ);
    Slocated=find(Mj4countS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    countPMj4=zeros(1,1);
    countTMj4=zeros(1,1);
    countP=zeros(1,1);
    countPbegin=zeros(1,1);
    countPend=zeros(1,1);
    countT=zeros(1,1);
    countTbegin = zeros(1,1);
    countTend = zeros(1,1);
    windowSize=100;
    %%%%%%%%%%%%%%%%%%%%%%% P����� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Rlocated Qlocated ���ڳ߶�4�µ�����
    for i=2:length(Rlocated)
        flag=0;
        mark4=0;
        RRinteral=Rlocated(i)-Rlocated(i-1);
        for j=1:5:(RRinteral*2/3)
            % windowEnd=Rlocated(i)-30-j;
            windowEnd=Qlocated(i)-j;
            windowBegin=windowEnd-windowSize;
            if windowBegin<Rlocated(i-1)+RRinteral/3
                break;
            end
            %���ڵļ���Сֵ
            %windowposi=Mj4.*(Mj4>0);
            %windowthposi=(max(Mj4(windowBegin:windowBegin+windowSize/2))+max(Mj4(windowBegin+windowSize/2+1:windowEnd)))/2;
            [windowMax,maxindex]=max(Mj4(windowBegin:windowEnd));
            [windowMin,minindex]=min(Mj4(windowBegin:windowEnd));
            if minindex < maxindex &&((maxindex-minindex)<windowSize*2/3)&&windowMax>0.01&&windowMin<-0.1
                flag=1;
                mark4=round((maxindex+minindex)/2+windowBegin);
                countPMj4(mark4)=1;
                countP(mark4-20)=1;
                countPbegin(windowBegin+minindex-20)=-1;
                countPend(windowBegin+maxindex-20)=-1;
            end
            if flag==1
                break;
            end
        end
        if mark4==0&&flag==0 %����û��P������R������1/3����������ֵ-1
            mark4=round(Rlocated(i)-RRinteral/3);
            countP(mark4-20)=-1;
        end
    end
    % figure(60000)
    %  plot(countPMj4,'g');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% T����� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear windowBegin windowEnd maxindex minindex windowMax windowMin mark4 RRinteral;
    
    windowSizeQ=100;
    for i=1:length(Rlocated)-1;
        flag=0;
        mark5=0;
        RRinteral=Rlocated(i+1)-Rlocated(i);
        for j=1:5:(RRinteral*2/3)
            % windowBegin=Rlocated(i)+30+j;
            windowBegin=Slocated(i)+j;
            windowEnd  =windowBegin+windowSizeQ;
            if windowEnd >Rlocated(i+1)-RRinteral/4
                break;
            end
            %%%%%�����ڵļ���Сֵ
            [windowMax,maxindex]=max(Mj4(windowBegin:windowEnd));
            [windowMin,minindex]=min(Mj4(windowBegin:windowEnd));
            if minindex < maxindex &&((maxindex-minindex)<windowSizeQ)&&windowMax>0.1&&windowMin<-0.1
                flag=1;
                mark5=round((maxindex+minindex)/2+windowBegin);
                countTMj4(mark5)=1;
                countT(mark5-20)=1;%�ҵ�T����ֵ��
                %%%%%ȷ��T����ʼ����յ�
                countTbegin(windowBegin+minindex-20)=-1;
                countTend(windowBegin+maxindex-20)=-1;
            end
            if flag==1
                break;
            end
        end
        if mark5==0 %����û��T������R���Ҽ��1/3����������ֵ-2
            mark5=round(Rlocated(i)+ RRinteral/3);
            countT(mark5)=-2;
        end
    end
    %plot(countTMj4,'g');
    %hold off;
    figure(9000);
    plot(ecgdata(0*points+1:1*points)),grid on,axis tight,axis([1,points,-2,5]);
    title('ECGѶ�ŵ�QRS������P�����');
    hold on
    plot(countR,'m'); %���R����ֵ
    plot(countQ,'r');
    plot(countS,'g');
    plot(countPbegin,'y');
    plot(countPend,'y');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4);
    plot(ecgdata(0*points+1:1*points)),grid on,axis tight,axis([1,points,-2,5]);
    title('ECGѶ�ŵĸ����������β��μ��');
    hold on
    plot(countR,'m'); %���R����ֵ
    plot(countQ,'r');
    plot(countS,'g');
    
    %%%%%%%%%%% ԭʼR��Ԥ������ж�ֵ�����룬���޷���ȷ��ʾ��R��λ�� %%%%%%%%%%%%%%
    % for i=1:Rnum
    %     if R_result(i)==0;
    %         break
    %     end
    %     plot(R_result(i),ecgdata(R_result(i)),'bo','MarkerSize',10,'MarkerEdgeColor','g');
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % plot(countP,'r');  %���P����ֵ
    % plot(countT,'r');  %���T����ֵ
    plot(countPbegin,'y');
    plot(countPend,'y');
    plot(countTbegin,'k');
    plot(countTend,'k');
    
    hold off
    
    %%
    %%%%%%%%%%%%%%%%%%%%%% �����Ѷ�� %%%%%%%%%%%%%%%%%%%%%%
    % figure(5);
    % plot(ecgdata),grid on,axis tight,axis([1,points,-2,5])
    
    %%%%%%%%%%%%%%%%%%%%%% �����Ѷ�� %%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%���С��ת�����Pwave��Rwave��Twave������λ����΢ƫ��%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%�ڽ�ʱ�����д���ԭʼѶ�ţ���Ѱ���ֵȷ����ֵλ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%P����λ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_peak_y=[];
    P_peak_x=[];
    y=ecgdata;
    i=1;
    num=1;
    Pbegin=-1;
    temp = 0;
    while (i<=length(y))
        if (i<=length(countPbegin))
            if(countPbegin(i)==-1)
                Pbegin=i;
                temp = i;
            end
        end
        if(Pbegin>-1)
            if (y(temp)<=y(i))
                temp = i;
            end
            if (countPend(i)==-1)
                P_peak_x(num)=temp;
                P_peak_y(num)=y(temp);
                num=num+1;
                Pbegin=-1;
            end
        end
        i=i+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% P����ͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(6);
    % plot(ecgdata),grid on,axis tight,axis([1,points,-2,5])
    % hold on
    % stem(P_peak_x,P_peak_y);
    % set(gcf,'units','normalized','position',[0,0.2,0.5,0.2]);
    % hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% R����λ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_peak_y=[];
    R_peak_x=[];
    y=ecgdata;
    i=1;
    num=1;
    Rbegin=1;
    temp = 0;
    while (i<=length(y))
        if (i<=length(countQ))
            if(countQ(i)==-1)
                Rbegin=i;
                temp = i;
            end
        end
        if(Rbegin>1)
            if (y(temp)<=y(i))
                temp = i;
            end
            if (countS(i)==-1)
                R_peak_x(num)=temp;
                R_peak_y(num)=y(temp);
                num=num+1;
                Rbegin=1;
            end
        end
        i=i+1;
    end
    % if(length(R_peak_x)~=length(countR))
    %     R_peak_x=countR;
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%% R����ͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(7);
    % plot(ecgdata),grid on,axis tight,axis([1,points,-2,5])
    % hold on
    % stem(R_peak_x,R_peak_y);
    % set(gcf,'units','normalized','position',[0,0.2,0.5,0.2]);
    % hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%% Q����λ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q_peak_y=[];
    Q_peak_x=[];
    y=ecgdata;
    i=1;
    num=1;
    Qbegin=1;
    temp = 0;
    while (i<=length(y))
        if (i<=length(countQ))
            if(countQ(i)==-1)
                Qbegin=i;
                temp = i;
            end
        end
        if(Qbegin>1)
            if (y(temp)>=y(i))
                temp = i;
            end
            if (countR(i)==1)
                Q_peak_x(num)=temp;
                Q_peak_y(num)=y(temp);
                num=num+1;
                Qbegin=1;
            end
        end
        i=i+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% Q����ͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(8);
    % plot(ecgdata),grid on,axis tight,axis([1,points,-2,5])
    % hold on
    % plot(Q_peak_x,Q_peak_y,'+');
    % set(gcf,'units','normalized','position',[0,0.2,0.5,0.2]);
    % hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%% S����λ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_peak_y=[];
    S_peak_x=[];
    y=ecgdata;
    i=1;
    num=1;
    Sbegin=1;
    temp = 0;
    while (i<=length(y))
        if (i<=length(countR))
            if(countR(i)==1)
                Sbegin=i;
                temp = i;
            end
        end
        if(Sbegin>1)
            if (y(temp)>=y(i))
                temp = i;
            end
            if (countS(i)==-1)
                S_peak_x(num)=temp;
                S_peak_y(num)=y(temp);
                num=num+1;
                Sbegin=1;
            end
        end
        i=i+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% S����ͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(9);
    % plot(ecgdata),grid on,axis tight,axis([1,points,-2,5])
    % hold on
    % plot(S_peak_x,S_peak_y,'^');
    % set(gcf,'units','normalized','position',[0,0.2,0.5,0.2]);
    % hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%T����λ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_peak_y=[];
    T_peak_x=[];
    y=ecgdata;
    i=1;
    num=1;
    Tbegin=-1;
    temp = 0;
    while (i<=length(y))
        if (i<=length(countTbegin))
            if(countTbegin(i)==-1)
                Tbegin=i;
                temp = i;
            end
        end
        if(Tbegin>-1)
            if (y(temp)<=y(i))
                temp = i;
            end
            if (countTend(i)==-1)
                T_peak_x(num)=temp;
                T_peak_y(num)=y(temp);
                num=num+1;
                Tbegin=-1;
            end
        end
        i=i+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% T����ͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(10);
    % plot(ecgdata),grid on,axis tight,axis([1,points,-2,5])
    % hold on
    % stem(T_peak_x,T_peak_y,'k');
    % set(gcf,'units','normalized','position',[0,0.2,0.5,0.2]);
    % hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PQRS��ͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(16+fix((te-1)/3));
    subplot(3,1,mod(te-1,3)+1)
    xvalue = linspace(0,points/sfreq,points);
    plot(xvalue,ecgdata),grid on,axis auto,ylabel('voltage(mV)'),xlabel('time(s)'),title(strcat(stringname," lead ",leadname(te)))
    hold on
    stem(P_peak_x/sfreq,P_peak_y,'o');
    plot(Q_peak_x/sfreq,Q_peak_y,'r+');
    stem(R_peak_x/sfreq,R_peak_y,'m');
    plot(S_peak_x/sfreq,S_peak_y,'g^');
    stem(T_peak_x/sfreq,T_peak_y,'k');
    hold off
    if ((mod(te-1,3)+1) == 3)
        saveas(16+fix((te-1)/3),strcat(stringname,'_',int2str(fix((te-1)/3))))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%% P-P interval �� R-R interval ��ֵ�Ƚ�%%%%%%%%%%%%%%%%%%%%
    %{
%%%%%%%%%%%%%%%%%%%%%%%%%%% P-P interval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPinter=[];
j=1;
for i=1:length(P_peak_x)-1
    PPinter(j)=((P_peak_x(i+1)-P_peak_x(i)))/sfreq;
    j=j+1;
end

HPV=mean(PPinter);           % PPintervalȡƽ��ֵ
 
samplerate1=(1/mean(PPinter))*60;
NN=length(PPinter);
% figure(100),plot(PPinter),title('Heart rate variability')
% xlabel('beats'),ylabel('P-P interval (s)')
figure(17),hist(PPinter,100),title(['HR = ' num2str(mean(PPinter)) '(s)' '  Heart rate = ' num2str(samplerate1,2) '(min)'])
xlabel('P-P interval (s)'),ylabel('count beats')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% R-R interval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RRinter=[];
j=1;
for i=1:length(R_peak_x)-1
    RRinter(j)=((R_peak_x(i+1)-R_peak_x(i)))/sfreq;
    j=j+1;
end

HRV=mean(RRinter);           % RRintervalȡƽ��ֵ��
 
samplerate2=(1/mean(RRinter))*60;
NN=length(RRinter);
% figure(100),plot(RRinter),title('Heart rate variability')
% xlabel('beats'),ylabel('R-R interval (s)')
figure(18),hist(RRinter,100),title(['HR = ' num2str(mean(RRinter)) '(s)' '  Heart rate = ' num2str(samplerate2,3) '(min)'])
xlabel('R-R interval (s)'),ylabel('count beats')
    %}
end
%%
%FUNCTIONS
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