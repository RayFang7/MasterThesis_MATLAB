%%
clc; clear;close all;
%%%%%%%%%%%讀取數據%%%%%%%%%%%%%%
%------ SPECIFY DATA ------------------------------------------------------
%%文件名稱
stringname='I01';
path='H:\Dataset\st-petersburg-incart-12-lead-arrhythmia-database-1.0.0\files\'; %原始訊號切割5s、圖片存檔位置
points=21600;
HEADERFILE= strcat(stringname,'.hea');     % header-file in text format
ATRFILE= strcat(stringname,'.atr');        % attributes-file in binary format
DATAFILE=strcat(stringname,'.dat');        % data-file
SAMPLES2READ=points;                       % 讀取數量 number of samples to be read
                                           % in case of more than one signal:
                                           % 2*SAMPLES2READ samples are read


%------ LOAD HEADER DATA --------------------------------------------------
fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE);
signalh= fullfile(PATH, HEADERFILE);
fid1=fopen(signalh,'r');
z= fgetl(fid1); 
A= sscanf(z, '%*s %d %d %d',[1,3]);
nosig= A(1);  % number of signals訊號數量
sfreq=A(2);   % sample rate of data取樣頻率
clear A;
for k=1:nosig
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
    dformat(k)= A(1);           % format; here only 212 is allowed
    gain(k)= A(2);              % number of integers per mV
    bitres(k)= A(3);            % bitresolution
    zerovalue(k)= A(4);         % integer value of ECG zero point
    firstvalue(k)= A(5);        % first integer value of signal (to test for errors)
end;
fclose(fid1);
clear A;

%------ LOAD BINARY DATA --------------------------------------------------
if dformat~= [212,212], error('this script does not apply binary formats different to 212.'); end; %不適用於二進制
signald= fullfile('C:','Users','AMD_107408104','Desktop','eddyhsu_paper','MIT_BIH_database', DATAFILE);            % data in format 212
fid2=fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit 
fclose(fid2);                                 % 矩陣有3行，每行8位，= 2 * 12位
M2H= bitshift(A(:,2), -4);
M1H= bitand(A(:,2), 15);
PRL=bitshift(bitand(A(:,2),8),9);     % sign-bit
PRR=bitshift(bitand(A(:,2),128),5);   % sign-bit
M( : , 1)= bitshift(M1H,8)+ A(:,1)-PRL;
M( : , 2)= bitshift(M2H,8)+ A(:,3)-PRR;
if M(1,:) ~= firstvalue, error('inconsistency in the first bit values'); end;  %第一位值不一致
switch nosig
case 2
    M( : , 1)= (M( : , 1)- zerovalue(1))/gain(1);
    M( : , 2)= (M( : , 2)- zerovalue(2))/gain(2);
    TIME=(0:(SAMPLES2READ-1))/sfreq;
case 1
    M( : , 1)= (M( : , 1)- zerovalue(1));
    M( : , 2)= (M( : , 2)- zerovalue(1));
    M=M';
    M(1)=[];
    sM=size(M);
    sM=sM(2)+1;
    M(sM)=0;
    M=M';
    M=M/gain(1);
    TIME=(0:2*(SAMPLES2READ)-1)/sfreq;
otherwise  % this case did not appear up to now!
    % here M has to be sorted!!!
    disp('Sorting algorithm for more than 2 signals not programmed yet!');
end;
clear A M1H M2H PRR PRL;
fprintf(1,'\\n$> LOADING DATA FINISHED \n');

%------ LOAD ATTRIBUTES DATA ----------------------------------------------
atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
fid3=fopen(atrd,'r');
A= fread(fid3, [2, inf], 'uint8')';
fclose(fid3);
ATRTIME=[];
ANNOT=[];
sa=size(A);
saa=sa(1);
i=1;
while i<=saa
    annoth=bitshift(A(i,2),-2);
    if annoth==59
        ANNOT=[ANNOT;bitshift(A(i+3,2),-2)];
        ATRTIME=[ATRTIME;A(i+2,1)+bitshift(A(i+2,2),8)+...
                bitshift(A(i+1,1),16)+bitshift(A(i+1,2),24)];
        i=i+3;
    elseif annoth==60
        % nothing to do!
    elseif annoth==61
        % nothing to do!
    elseif annoth==62
        % nothing to do!
    elseif annoth==63
        hilfe=bitshift(bitand(A(i,2),3),8)+A(i,1);
        hilfe=hilfe+mod(hilfe,2);
        i=i+hilfe/2;
    else
        ATRTIME=[ATRTIME;bitshift(bitand(A(i,2),3),8)+A(i,1)];
        ANNOT=[ANNOT;bitshift(A(i,2),-2)];
   end;
   i=i+1;
end;
ANNOT(length(ANNOT))=[];       % last line = EOF (=0)
ATRTIME(length(ATRTIME))=[];   % last line = EOF
clear A;
ATRTIME= (cumsum(ATRTIME))/sfreq;
ind= find(ATRTIME <= TIME(end));
ATRTIMED= ATRTIME(ind);
ANNOT=round(ANNOT);
ANNOTD= ANNOT(ind);

%------ DISPLAY DATA ------------------------------------------------------
% y=detrend(M(:,1));
% plot(TIME, y ,'b');
% 
% 
% figure(1); clf, box on, hold on ;grid on ;
% plot(TIME, M(:,1),'b');
% if nosig==2
%     figure(2); clf, box on, hold on ;grid on ;
%     plot(TIME, M(:,2),'b');
% end;
% for k=1:length(ATRTIMED)
%     text(ATRTIMED(k),0,num2str(ANNOTD(k)));
% end;
% xlim([TIME(1), TIME(end)]);
% xlabel('Time / s'); ylabel('Voltage / mV');
% string=['ECG signal ',DATAFILE];
% title(string);



fprintf(1,'\\n$> DISPLAYING DATA FINISHED \n');
% -------------------------------------------------------------------------
fprintf(1,'\\n$> ALL FINISHED \n');
% hold off

%%%%%%%%%%%%%%%%%%%去除雜訊噪聲和基線漂移%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%小波轉換%%%
level=8; wavename='coif5';
ecgdata=M(:,1);

figure(3)
plot(ecgdata(1:points));grid off ;axis tight;
axis([1,points,-4,4]);
title('原始ECG訊號');
xlabel('資料點位置');
ylabel('電壓(mV)')



%%%%%%%%%%進行小波轉換8層
[C,L]=wavedec(ecgdata,level,wavename);
%%%%%%%提取小波尺度系數，
A1=appcoef(C,L,wavename,1);
A2=appcoef(C,L,wavename,2);
A3=appcoef(C,L,wavename,3);
A4=appcoef(C,L,wavename,4);
A5=appcoef(C,L,wavename,5);
A6=appcoef(C,L,wavename,6);
A7=appcoef(C,L,wavename,7);
A8=appcoef(C,L,wavename,8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(100),plot(A1)
% figure(101),plot(A2)
% figure(102),plot(A3)
% figure(103),plot(A4)
% figure(104),plot(A5)
% figure(105),plot(A6)
% figure(106),plot(A7)
% figure(107),plot(A8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%提取細節系數
D1=detcoef(C,L,1);
D2=detcoef(C,L,2);
D3=detcoef(C,L,3);
D4=detcoef(C,L,4);
D5=detcoef(C,L,5);
D6=detcoef(C,L,6);
D7=detcoef(C,L,7);
D8=detcoef(C,L,8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(108),plot(D1)
% figure(109),plot(D2)
% figure(110),plot(D3)
% figure(111),plot(D4)
% figure(112),plot(D5)
% figure(113),plot(D6)
% figure(114),plot(D7)
% figure(115),plot(D8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%重構
A8=zeros(length(A8),1); %去除基線漂移,8層低頻信息
RA7=idwt(A8,D8,wavename);
RA6=idwt(RA7(1:length(D7)),D7,wavename);
RA5=idwt(RA6(1:length(D6)),D6,wavename);
RA4=idwt(RA5(1:length(D5)),D5,wavename);
RA3=idwt(RA4(1:length(D4)),D4,wavename);

RA2=idwt(RA3(1:length(D3)),D3,wavename);
D2=zeros(length(D2),1); %去除高頻噪聲，2層高頻噪聲
RA1=idwt(RA2(1:length(D2)),D2,wavename);
D1=zeros(length(D1),1);%去除高頻噪聲，1層高頻噪聲
DenoisingSignal=idwt(RA1,D1,wavename);

DenoisingSignal_xx=0;
for i=1:length(DenoisingSignal)
    DenoisingSignal_xx=DenoisingSignal_xx+1;
    DenoisingSignal_x(i)=DenoisingSignal_xx;
end
DenoisingSignal_x=DenoisingSignal_x';

figure(4);
plot(DenoisingSignal);
title('去除高頻訊號與基線飄移的ECG訊號'); grid on; axis tight;axis([1,points,-4,4]);


%---標記R波
% R_xpoint=round(ATRTIMED*360+1);
% R_xpoint=R_xpoint(2:length(R_xpoint));
% for a=1:length(R_xpoint)
%   R_ypoint(a,1)=DenoisingSignal(R_xpoint(a,1),1);
% end
% figure(10)
% plot(DenoisingSignal);
% title('R點標記');
% grid on; axis tight;axis([1,points,-4,4]);
% hold on;
% plot(R_xpoint,R_ypoint,'o');
%% 原始訊號切割5s、圖片存檔位置

% for j=1:360
%     if j==1
%     a1(j,:)=[1,1800];
%     elseif j>1
%     a1(j,:)=[a1(j-1,2)+1,a1(j-1,2)+1800];    
%     end
% end
% 
% for i=1:360
%     figure(5)
%     plot(ecgdata(1:points),'k') ;
%     axis([a1(i,:),-4,4]);
%     axis off;
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 15])
%     print('-dpng',[path ,'5s_',num2str(i),'.png'],'-r200');
% end



%% 原始訊號，用harris corner角點分析找所有端點



[corners_x,corners_y,I2]=harriscornerdetect(path2);


figure(6)
imshow(I2);hold on;
title('Harris所有角點標記');

plot(corners_x,corners_y,'rx'); %找所有端點
legend('harris角點');

% figure(106)
% imhist(I2);


%% R波位置標記
% flag=0;
% minn=min(corners.Location(:,2));
% maxx=max(corners.Location(:,2));
% average=(minn+maxx)*1.3/3;
% for i=1:length(corners.Location(:,2))
%     if corners.Location(i,2)<average             
%         flag=1; 
%         loca(i,1)=flag;             
%     elseif corners.Location(i,2)>=average
%         flag=0;
%         loca(i,1)=flag;
%     end 
% end
% for k=1:length(loca)
%     if k+1<=length(loca)
%         if loca(k+1,1)==1 && loca(k+1,1)==loca(k,1)
%         loca(k,1)=1;
%         loca(k+1,1)=0;
%         
%         end   
%     end
% end
% loca2=find(loca==1);
% for j=1:length(loca2)
%     R(j,:)=corners.Location(loca2(j),:);
% end
% % figure(7)
% % imshow(I2);hold on;
% % title('R波標記');
% % plot(R(:,1),R(:,2),'o'); %R波位置標記
% % legend('R波');
%% R波位置標記(test1)

[corners_x2,corners_y2]=localmin(corners_x,corners_y); % 找y軸局部最小值法(找R波)

figure(1111)
imshow(I2);hold on;
title('Harris所有角點標記(y軸局部最小值法)');
plot(corners_x2,corners_y2,'rx'); %找所有端點
legend('harris角點');



threshold1=findthreshold(corners_y2);  % 直方圖(迭代法選擇閾值)

% figure('Name','corners_y2')
% aa1=bar(corners_y2);

    a2=find(corners_y2<threshold1); %找出小於閾值的點
     for i=1:length(a2)
        corners_x3(i,1)=corners_x2(a2(i));
        corners_y3(i,1)=corners_y2(a2(i));
     end
% figure(77)
% imshow(I2);hold on;
% title('R波標記');
% plot(corners_x3(:,1),corners_y3(:,1),'o'); %R波位置標記
% legend('R波');
       
% figure('Name','corners_y3')
% aa2=bar(corners_y3);

if length(a2)>8    %設定R波數量在8以下
    threshold1=findthreshold(corners_y3);
    clear a2;
    a2=find(corners_y3<threshold1);
    for j=1:length(a2)
       corners_x4(j,1)=corners_x3(a2(j));
       corners_y4(j,1)=corners_y3(a2(j)); 
    end     
       if length(a2)>8
          threshold1=findthreshold(corners_y4);
          clear a2;
          a2=find(corners_y4<threshold1);
          for k=1:length(a2)
          corners_x5(k,1)=corners_x4(a2(k));
          corners_y5(k,1)=corners_y4(a2(k)); 
          end 
          R(:,2)=corners_x5;
          R(:,3)=corners_y5;

       elseif length(a2)<=8
          R(:,2)=corners_x4;
          R(:,3)=corners_y4;
       end
     
elseif length(a2)<=8
    R(:,2)=corners_x3;
    R(:,3)=corners_y3;
end


% [AA1,AA2]=kmeans(R(:,2),5);
% for k=1:length(AA2)
%     for m=1:length(R(:,2))        
%       newA(m,k)=AA2(k,1)-R(m,2);     
%     end
% end
% newA=abs(newA);
% [aa1,aa2]=min(newA);
% aa2=aa2';
% aa3=sort(aa2,'ascend');
% for n=1:length(aa3)
% aa4(n,2:3)=R(aa3(n,1),2:3);
% end




% figure('Name','corners_y4')
% aa3=bar(corners_y4);
% 
% figure('Name','corners_y5')
% aa4=bar(corners_y5);

figure(7)
imshow(I2);hold on;
title('R波標記');
plot(R(:,2),R(:,3),'o'); %R波位置標記
legend('R波');



%% Q波位置標記
% for i=1:length(loca2)
%     if loca2(i)==1
%      matrix5(i,1)=0;
%     elseif loca2(i)==2
%      matrix5(i,1)=0;
%     elseif loca2(i)==3
%      matrix5(i,1)=0;
%     elseif loca2(i)>3
%      matrix5(i,1)=max(corners.Location(loca2(i)-Q_forward_parameter1:loca2(i)-Q_forward_parameter2,2));
%     end
% end
% for l=1:length(matrix5)
%     if matrix5(l,1)==0
%        matrix5(l,1)=0;
%     elseif matrix5(l,1)==matrix5(l-1,1)
%        matrix5(l,1)=0; 
%     elseif matrix5(l,1)~=0
%        matrix5(l,1)=matrix5(l,1);
%     end
% end
% corners_y=corners.Location(:,2);
% for j=1:length(matrix5)
%     if matrix5(j,1)==0
%     matrix6(j,1)=0;
%     elseif matrix5(j,1)~=0
%     a=[];
%     a=find(corners_y==matrix5(j,1));    
%     matrix6(j,1)=a(1,1);
%     end
% end
% for k=1:length(matrix6)
%     if matrix6(k,1)==0
%     Q(k,:)=[0,0];
%     elseif matrix6(k,1)~=0
%     Q(k,:)=corners.Location(matrix6(k,1),:);
%     end
% end
% figure(8)
% imshow(I2);hold on;
% title('Q波標記');
% plot(Q(:,1),Q(:,2),'g*'); %Q波位置標記
% legend('Q波');
%% Q波位置標記(test1)
for i=1:length(R(:,1))
   R(i,1)=find(corners_x==R(i,2));
end

for j=1:length(R(:,1))
    if R(j,1)==1
        Q(j,1)=0;
        Q(j,3)=0;
    elseif R(j,1)==2
        Q(j,1)=0;
        Q(j,3)=0;
    elseif R(j,1)==3
        Q(j,1)=0;
        Q(j,3)=0;
    elseif R(j,1)>3
        [b1(j),b2(j)]=max(corners_y(R(j,1)-Q_forward_parameter1:R(j,1)-Q_forward_parameter2)); %Q_forward_parameter1=3;Q_forward_parameter2=1
        if b2(j)==3
        Q(j,1)=R(j,1)-1;
        Q(j,3)=b1(j);
        elseif b2(j)==2
        Q(j,1)=R(j,1)-2;    
        Q(j,3)=b1(j);
        elseif b2(j)==1
        Q(j,1)=R(j,1)-3; 
        Q(j,3)=b1(j);
        end
    end
end
for k=1:length(Q(:,1))
    if Q(k,1)~=0
       Q(k,2)=corners_x(Q(k,1));
    elseif Q(k,1)==0
       Q(k,2)=0;  
    end
end


% figure(8)
% imshow(I2);hold on;
% title('Q波標記');
% plot(Q(:,2),Q(:,3),'g*'); %Q波位置標記
% legend('Q波');


%% S波位置標記
% for i=1:length(loca2)
%      if loca2(i)+S_backward_parameter2<=length(corners)
%      matrix7(i,1)=max(corners.Location(loca2(i)+S_backward_parameter1:loca2(i)+S_backward_parameter2,2));
%      elseif loca2(i)+S_backward_parameter2>length(corners)
%      matrix7(i,1)=0;
%      end
% end
% 
% corners_y=corners.Location(:,2);
% for j=1:length(matrix7)
%     if matrix7(j,1)==0
%     matrix8(j,1)=0;
%     elseif matrix7(j,1)~=0
%     b=[];
%     b=find(corners_y==matrix7(j,1));
%     matrix8(j,1)=b(1,1);
%     end
% end
% for k=1:length(matrix8)
%     if matrix8(k,1)==0
%     S(k,:)=[0,0];
%     elseif matrix8(k,1)~=0
%     S(k,:)=corners.Location(matrix8(k,1),:);
%     end
% end
% % figure(9)
% % imshow(I2);hold on;
% % title('S波標記');
% % plot(S(:,1),S(:,2),'bo'); %S波位置標記
% % legend('S波');

%% S波位置標記(test1)
for i=1:length(R(:,1))                                          %S_backward_parameter1=1;  S_backward_parameter2=3
     if R(i,1)+S_backward_parameter2<=length(corners_x)
        [c1(i),c2(i)]=max(corners_y(R(i,1)+S_backward_parameter1:R(i,1)+S_backward_parameter2));
        if c2(i)==1
           S(i,1)=R(i,1)+1;
           S(i,3)=c1(i);
        elseif c2(i)==2
           S(i,1)=R(i,1)+2;
           S(i,3)=c1(i);
        elseif c2(i)==3
           S(i,1)=R(i,1)+3;
           S(i,3)=c1(i);
        end
     elseif R(i,1)+S_backward_parameter2>length(corners_x)
     S(i,1)=0;
     S(i,3)=0;
     end
end
for j=1:length(S(:,1))
    if S(j,1)~=0
    S(j,2)=corners_x(S(j,1));
    elseif S(j,1)==0
    S(j,2)=0;  
    end
end


% figure(9)
% imshow(I2);hold on;
% title('S波標記');
% plot(S(:,2),S(:,3),'bo'); %S波位置標記
% legend('S波');



%% S波結束點標記
% for i=1:length(matrix8)
%    
%      matrix15(i,1)=min(corners.Location(matrix8(i)+1:matrix8(i)+2,2));
%     
% end
% corners_y=corners.Location(:,2);
% for j=1:length(matrix15)
%     if matrix15(j,1)==0
%     matrix16(j,1)=0;
%     elseif matrix15(j,1)~=0
%     matrix16(j,1)=find(corners_y==matrix15(j,1));
%     end
% end
% for k=1:length(matrix16)
%     if matrix16(k,1)==0
%     S_end(k,:)=[0,0];
%     elseif matrix16(k,1)~=0
%     S_end(k,:)=corners.Location(matrix16(k,1),:);
%     end
% end
% figure(10)
% imshow(I2);hold on;
% title('S波結束點標記');
% plot(S_end(:,1),S_end(:,2),'m+'); %S波結束點位置標記
% legend('S波結束點');
%% Q波起始點位置標記
% for i=1:length(matrix6)
%      if matrix6(i)<=6
%      matrix11(i,1)=0;
%      elseif matrix6(i)>6  
%      matrix11(i,1)=min(corners.Location(matrix6(i)-Qstart_forward_parameter1:matrix6(i)-Qstart_forward_parameter2,2));
%      end
% end
% corners_y=corners.Location(:,2);
% for j=1:length(matrix11)
%     if matrix11(j,1)==0
%     matrix12(j,1)=0;
%     elseif matrix11(j,1)~=0
%     matrix12(j,1)=find(corners_y==matrix11(j,1));
%     end
% end
% for k=1:length(matrix12)
%     if matrix12(k,1)==0
%     Q_start(k,:)=[0,0];
%     elseif matrix12(k,1)~=0
%     Q_start(k,:)=corners.Location(matrix12(k,1),:);
%     end
% end
% % figure(11)
% % imshow(I2);hold on;
% % title('Q波起始點標記');
% % plot( Q_start(:,1), Q_start(:,2),'r+'); %Q波起始點位置標記
% % legend('Q波起始點');

%% T波起始點位置標記
% for i=1:length(matrix10)
%     if matrix10(i)~=0
%     matrix13(i,1)=max(corners.Location(matrix10(i)-Tstart_forward_parameter1:matrix10(i)-Tstart_forward_parameter2,2)); 
%     else
%     matrix13(i,1)=0;    
%     end
% end
% corners_y=corners.Location(:,2);
% for j=1:length(matrix13)
%     if matrix13(j,1)==0
%     matrix14(j,1)=0;
%     elseif matrix13(j,1)~=0
%     c=[];
%     c=find(corners_y==matrix13(j,1));
%     matrix14(j,1)=c(1,1);
%     end
% end
% for k=1:length(matrix14)
%     if matrix14(k,1)==0
%     T_start(k,:)=[0,0];
%     elseif matrix14(k,1)~=0
%     T_start(k,:)=corners.Location(matrix14(k,1),:);
%     end
% end
% figure(13)
% imshow(I2);hold on;
% title('T波起始點標記');
% plot(T_start(:,1),T_start(:,2),'ro'); %T波起始點位置標記
% legend('T波起始點');
%% P波位置標記
%  for j=1:length(matrix6)
%      if matrix6(j)==0
%      matrix(j,1)=0;
%      elseif matrix6(j)~=0
%         if matrix6(j)-P_forward_parameter1>0
%         matrix(j,1)=min(corners.Location(matrix6(j)-P_forward_parameter1:matrix6(j)-P_forward_parameter2,2));
%         elseif matrix6(j)-P_forward_parameter1<=0
%         matrix(j,1)=corners.Location(1,2);
%         end
%      end
%  end
%  corners_y=corners.Location(:,2);
% for k=1:length(matrix)
%     if matrix(k,1)==0
%     matrix2(k,1)=0;
%     elseif matrix(k,1)~=0
%     matrix2(k,1)=find(corners_y==matrix(k,1));
%     end
% end
% for l=1:length(matrix2)
%     if matrix2(l,1)==0
%     P(l,:)=[0,0];
%     elseif matrix2(l,1)~=0
%     P(l,:)=corners.Location(matrix2(l,1),:);
%     end
% end
% % figure(14)
% % imshow(I2);hold on;
% % title('P波標記');
% % plot(P(:,1),P(:,2),'rx'); %P波位置標記
% % legend('P波');
%% P波位置標記(test1)
%從corners_y2局部y軸最小值找P波
%從R波往前找

for i=1:length(R(:,1))                     %找R波在corners_x2的位置dd1
    dd1(i,1)=find(corners_x2==R(i,2));
end
for j=1:length(dd1(:,1))
    if dd1(j,1)-7>0
    [e1(j),e2(j)]=min(corners_y2(dd1(j,1)-7:dd1(j,1)-1));
    P(j,1)=dd1(j,1)-8+e2(j);
    P(j,3)=e1(j);
    elseif dd1(j,1)-7==0 || dd1(j,1)-7<0
        if dd1(j,1)-1~=0
          [e1(j),e2(j)]=min(corners_y2(1:dd1(j,1)-1));    
          P(j,1)=dd1(j,1)-dd1(j,1)+e2(j);
          P(j,3)=e1(j);
        elseif dd1(j,1)-1==0
          P(j,1)=0;
          P(j,3)=0;
        end
    end
end
for k=1:length(P(:,1))
    if P(k,1)~=0
    P(k,2)=corners_x2(P(k,1),1);
    elseif P(k,1)==0
    P(k,2)=0;    
    end
end

% figure(14)
% imshow(I2);hold on;
% title('P波標記');
% plot(P(:,2),P(:,3),'rx'); %P波位置標記
% legend('P波');
%% P波起始點標記
% for i=1:length(matrix2)
%     if matrix2(i)~=0
%         if matrix2(i)-Pstart_forward_parameter1>0
%         matrix3(i,1)=max(corners.Location(matrix2(i)-Pstart_forward_parameter1:matrix2(i)-Pstart_forward_parameter2,2));
%         elseif matrix2(i)-Pstart_forward_parameter1<=0
%         matrix3(i,1)=corners.Location(1,2);
%         end
%     elseif matrix2(i)==0
%      matrix3(i,1)=0;
%     end
% end
% corners_y=corners.Location(:,2);
% for j=1:length(matrix3)
%     if matrix3(j,1)==0
%     matrix4(j,1)=0;
%     elseif matrix3(j,1)~=0
%     matrix4(j,1)=find(corners_y==matrix3(j,1));
%     end
% end
% for k=1:length(matrix4)
%     if matrix4(k,1)==0
%     P_start(k,:)=[0,0];
%     elseif matrix4(k,1)~=0
%     P_start(k,:)=corners.Location(matrix4(k,1),:);
%     end
% end
% % figure(15)
% % imshow(I2);hold on;
% % title('P波起始點標記');
% % plot(P_start(:,1),P_start(:,2),'m+'); %P波起始點位置標記
% % legend('P波起始點');

%% P波起始點標記(test1)
%從corners_y2局部y軸最大值找P波起始點
%從P波往前找

for i=1:length(P(:,1))                     
    if P(i,1)-6>0
    [f1(i),f2(i)]=max(corners_y2(P(i,1)-6:P(i,1)-1));
    P_start(i,1)=P(i,1)-7+f2(i);
    P_start(i,3)=f1(i);
    elseif P(i,1)-6==0 || P(i,1)-6<0
        if P(i,1)-1>0
        [f1(i),f2(i)]=max(corners_y2(1:P(i,1)-1));
        P_start(i,1)=P(i,1)-P(i,1)+f2(i);
        P_start(i,3)=f1(i);
        elseif P(i,1)-1<=0
        P_start(i,1)=P(i,1);
        P_start(i,3)=P(i,3);
        end
    end
end
for k=1:length(P_start(:,1))
    if P_start(k,1)~=0
       P_start(k,2)=corners_x2(P_start(k,1),1);
    elseif P_start(k,1)==0
       P_start(k,2)=0;
    end
end

% figure(15)
% imshow(I2);hold on;
% title('P波起始點標記');
% plot(P_start(:,2),P_start(:,3),'m+'); %P波起始點位置標記
% legend('P波起始點');
%% T波位置標記
% for i=1:length(matrix8)
%    if matrix8(i)+T_backward_parameter2<=length(corners.Location)
%    matrix9(i,1)=min(corners.Location(matrix8(i)+T_backward_parameter1:matrix8(i)+T_backward_parameter2,2));
%    else
%    matrix9(i,1)=0;    
%    end
% end
% corners_y=corners.Location(:,2);
% for j=1:length(matrix9)
%     if matrix9(j,1)==0
%     matrix10(j,1)=0;
%     elseif matrix9(j,1)~=0
%     matrix10(j,1)=find(corners_y==matrix9(j,1));
%     end
% end
% for k=1:length(matrix10)
%     if matrix10(k,1)==0
%     T(k,:)=[0,0];
%     elseif matrix10(k,1)~=0
%     T(k,:)=corners.Location(matrix10(k,1),:);
%     end
% end
% % figure(12)
% % imshow(I2);hold on;
% % title('T波標記');
% % plot(T(:,1),T(:,2),'r^'); %T波位置標記
% % legend('T波');
%% T波位置標記(test1)

% %從corners_y2局部y軸最小值找T波
% %從R波往後找
% 
% for i=1:length(R(:,1))                     %找R波在corners_x2的位置dd1
%     dd1(i,1)=find(corners_x2==R(i,2));
% end
% 
% for j=1:length(dd1(:,1))
%     if dd1(j,1)+10<=length(corners_x2)
%     [d1(j),d2(j)]=min(corners_y2(dd1(j,1)+1:dd1(j,1)+10));
%     T(j,1)=dd1(j,1)+d2(j);
%     T(j,3)=d1(j);
%     elseif dd1(j,1)+10>length(corners_x2)
%     [d1(j),d2(j)]=min(corners_y2(dd1(j,1)+1:length(corners_x2)));    
%     T(j,1)=dd1(j,1)+d2(j);
%     T(j,3)=d1(j);
%     end
% end
% 
% 
% for k=1:length(T(:,1))
%     T(k,2)=corners_x2(T(k,1),1);
% end
% 
% % figure(12)
% % imshow(I2);hold on;
% % title('T波標記');
% % plot(T(:,2),T(:,3),'r^'); %T波位置標記
% % legend('T波');
%% T波位置標記(test2)

%從R波與下一個P波之間y軸最小值找T波


for i=1:length(R(:,1))                     %找R波在corners_x2的位置dd1
    dd1(i,1)=find(corners_x2==R(i,2));
end

for j=1:length(dd1(:,1))
    if j+1<=length(dd1(:,1))
    [d1(j),d2(j)]=min(corners_y2(dd1(j,1)+1:P_start(j+1,1)-1));
    T(j,1)=dd1(j,1)+d2(j);
    T(j,3)=d1(j);
    elseif j+1>length(dd1(:,1))     
        if dd1(j,1)+20<=length(corners_x2)
        [d1(j),d2(j)]=min(corners_y2(dd1(j,1)+1:dd1(j,1)+20)); 
        T(j,1)=dd1(j,1)+d2(j);
        T(j,3)=d1(j);
        elseif dd1(j,1)+20>length(corners_x2) && dd1(j,1)~=length(corners_x2)
        [d1(j),d2(j)]=min(corners_y2(dd1(j,1)+1:length(corners_x2)));
        T(j,1)=dd1(j,1)+d2(j);
        T(j,3)=d1(j);
        elseif dd1(j,1)==length(corners_x2)   
        T(j,1)=0;
        T(j,3)=0;    
        end 
    end
end


for k=1:length(T(:,1))
    if T(k,1)~=0
    T(k,2)=corners_x2(T(k,1),1);
    elseif T(k,1)==0
    T(k,2)=0;    
    end
end

% figure(12)
% imshow(I2);hold on;
% title('T波標記');
% plot(T(:,2),T(:,3),'r^'); %T波位置標記
% legend('T波');
%% 顯示全部標記位置
figure(16)
imshow(I2);hold on;
title('主要波形標記');
%P起始點、PQRST標記
plot(P_start(:,2),P_start(:,3),'m+',P(:,2),P(:,3),'bx',Q(:,2),Q(:,3),'g*',R(:,2),R(:,3),'o',S(:,2),S(:,3),'bo',T(:,2),T(:,3),'r^'); %主要波形標記
legend('P波起始點','P波','Q波','R波','S波','T波');
%QRS標記
% plot(Q(:,2),Q(:,3),'g*',R(:,2),R(:,3),'o',S(:,2),S(:,3),'bo'); %主要波形標記
% legend('Q波','R波','S波');
%P波起始點、PQRS標記
% plot(P_start(:,2),P_start(:,3),'m+',P(:,2),P(:,3),'bx',Q(:,2),Q(:,3),'g*',R(:,2),R(:,3),'o',S(:,2),S(:,3),'bo'); %主要波形標記
% legend('P波起始點','P波','Q波','R波','S波');
%% 自動切圖片-onewave
% P_x=round(P_start(:,2));
% y_top=min(R(:,3));
% y_top=y_top-100;
% y_1=max(Q(:,3));
% y_2=max(S(:,3));
% y_bottom=max(y_1,y_2);
% y_bottom=y_bottom+100;
% clear y_1 y_2;
% y_top=round(y_top);
% y_bottom=round(y_bottom);
% 
% for i=1:length(P_x)
%     if i>=1 && i<length(P_x)
%        if P_x(i)~=0
%          I3=I2(y_top:y_bottom,P_x(i):P_x(i+1));
%          figure(17)
%          imshow(I3);
%          print('-dpng',[path3 ,num2str(i),'.png'],'-r600');
%        elseif P_x(i)==0
%          I3=I2(y_top:y_bottom,P_x(i+1):P_x(i+2));
%          figure(17)
%          imshow(I3);
%          print('-dpng',[path3 ,num2str(i),'.png'],'-r600');
%        end
%     elseif i==length(P_x)
%          I3=I2(y_top:y_bottom,P_x(i):length(I2));
%          figure(17)
%          imshow(I3);
%          print('-dpng',[path3 ,num2str(i),'.png'],'-r600');
% 
%     end
% end
%% 自動切圖片-twowave
% P_x=round(P_start(:,1));
% y_top=min(R(:,2));
% y_top=y_top-100;
% y_1=max(Q(:,2));
% y_2=max(S(:,2));
% y_bottom=max(y_1,y_2);
% y_bottom=y_bottom+100;
% clear y_1 y_2;
% y_top=round(y_top);
% y_bottom=round(y_bottom);
% for i=1:length(P_x)
%     if i>=1 && i<length(P_x)
%        if P_x(i)~=0
%          I3=I2(y_top:y_bottom,P_x(i):P_x(i+2));
%          figure(17)
%          imshow(I3);
%          print('-dpng',[path4 ,num2str(i),'.png'],'-r600');
%        elseif P_x(i)==0
%          I3=I2(y_top:y_bottom,P_x(i+1):P_x(i+3));
%          figure(17)
%          imshow(I3);
%          print('-dpng',[path4 ,num2str(i),'.png'],'-r600');
%        end
%     elseif i==length(P_x)
%        break;
%     end
% end
%% 自動切圖片-threewave
% P_x=round(P_start(:,1));
% y_top=min(R(:,2));
% y_top=y_top-100;
% y_1=max(Q(:,2));
% y_2=max(S(:,2));
% y_bottom=max(y_1,y_2);
% y_bottom=y_bottom+100;
% clear y_1 y_2;
% y_top=round(y_top);
% y_bottom=round(y_bottom);
% for i=1:length(P_x)
%     if i>=1 && i<length(P_x)
%        if P_x(i)~=0
%          I3=I2(y_top:y_bottom,P_x(i):P_x(i+3));
%          figure(17)
%          imshow(I3);
%          print('-dpng',[path5 ,num2str(i),'.png'],'-r600');
%        elseif P_x(i)==0
%          I3=I2(y_top:y_bottom,P_x(i+1):P_x(i+4));
%          figure(17)
%          imshow(I3);
%          print('-dpng',[path5 ,num2str(i),'.png'],'-r600');
%        end
%     elseif i==length(P_x)
%        break;
%     end
% end
%% 手動切圖片(多張)
% P_x=[611;1032;1370;1867;2307;2697];
% y_top=min(R(:,2));
% y_top=y_top-100;
% y_1=max(Q(:,2));
% y_2=max(S(:,2));
% y_bottom=max(y_1,y_2);
% y_bottom=y_bottom+100;
% clear y_1 y_2;
% y_top=round(y_top);
% y_bottom=round(y_bottom);
% for i=1:length(P_x)
%     if i>=1 && i<length(P_x)
%        if P_x(i)~=0
%          I3=I2(y_top:y_bottom,P_x(i):P_x(i+1));
%          figure(17)
%          imshow(I3);hold on;
%          print('-dpng',[path ,num2str(i),'.png'],'-r600');
%        elseif P_x(i)==0
%          I3=I2(y_top:y_bottom,P_x(i+1):P_x(i+2));
%          figure(17)
%          imshow(I3);hold on;
%          print('-dpng',[path ,num2str(i),'.png'],'-r600');
%        end
%     elseif i==length(P_x)
%        break;
%     end
% end
%% 手動切圖片(單張)
% P_x=[1066;1228];
% y_top=min(R(:,2));
% y_top=y_top-100;
% y_1=max(Q(:,2));
% y_2=max(S(:,2));
% y_bottom=max(y_1,y_2);
% y_bottom=y_bottom+100;
% clear y_1 y_2;
% y_top=round(y_top);
% y_bottom=round(y_bottom);
% 
% I3=I2(y_top:y_bottom,P_x(1):P_x(2));
% figure(18)
% imshow(I3);hold on;
% print('-dpng',[path3 ,'1.png'],'-r600');
%% 手動標記(單張)
% figure(17)
% imshow(I2);hold on;
% title('主要波形標記');
% plot(P_start(:,2),P_start(:,3),'m+',P(:,2),P(:,3),'bx',Q(:,2),Q(:,3),'g*',R(:,2),R(:,3),'o',S(:,2),S(:,3),'bo',T(:,2),T(:,3),'r^',[1222.2886,1825],[1487.7445,1547],'r*',[1467.8813,2044],[1472.7618,1513],'b*'); %主要波形標記
% legend('P波起始點','P波','Q波','R波','S波','T波','VT起點','VT終點');
%%
%----------切割圖片--------------------%
% 分離點標記(一個波形)

% clear splitpoint_x splitpoint_y
% wave_atr=ANNOTD(2:end,1);
% for k1=1:length(R_xpoint)
%     if k1<length(R_xpoint)
%        splitpoint(k1,1)=round(abs(R_xpoint(k1,1)-R_xpoint(k1+1,1))/2);
%     else
%        break;
%     end
% end
% for k2=1:length(splitpoint)
%       splitpoint2(k2,1)=splitpoint(k2,1)+R_xpoint(k2,1);
% end
% splitpoint_x=splitpoint2;
% clear splitpoint2 splitpoint ;
% 
% for k3=1:length(splitpoint_x)
%     splitpoint_y(k3,1)=DenoisingSignal(splitpoint_x(k3,1),1);
% end
% 
% figure(6)
% plot(DenoisingSignal);
% title('分離點標記');
% grid on; axis tight;axis([1,points,-4,4]);hold on;
% plot(splitpoint_x,splitpoint_y,'r+');

%%
%----------切割圖片--------------------%
% 分離點標記(二個波形)
% splitpoint_x1=splitpoint_x(2:2:end);
% splitpoint_y1=splitpoint_y(2:2:end);
% figure(66)
% plot(DenoisingSignal);
% title('偶數分離點標記');
% grid on; axis tight;axis([1,points,-4,4]);hold on;
% plot(splitpoint_x1,splitpoint_y1,'r+');
% 
% 
% fprintf(1,'\\n$> DIVIDE PICTURES FINISHED \n');

%%
%---存圖片------------------------------------%
%----圖片大小(875x656)------------------------%
%-----切割1個波形----------------------------%
% wave_atr=ANNOTD(2:end,1);
% 
% for j=1:length(wave_atr)
%     if j==1
%        figure(7)
%        plot(DenoisingSignal,'k');
%        axis([1,splitpoint_x(j),-4,4]);
%        axis off
%        print('-dpng',[location1 ,num2str(j),'.png'],'-r150');
%     elseif j>1&&j<length(wave_atr)
%        figure(7)
%        plot(DenoisingSignal,'k');
%        axis([splitpoint_x(j-1),splitpoint_x(j),-4,4]);
%        axis off
%        print('-dpng',[location1 ,num2str(j),'.png'],'-r150');
%     elseif j==length(wave_atr)
%        figure(7)
%        plot(DenoisingSignal,'k');
%        axis([splitpoint_x(j-1),points,-4,4]);
%        axis off
%        print('-dpng',[location1 ,num2str(j),'.png'],'-r150');
%         
%     end
% end
% 
% fprintf(1,'\\n$> SAVE PICTURES FINISHED \n');
%%
%---存圖片------------------------------------%
%----圖片大小(875x656)------------------------%
%-----切割2個波形----------------------------%
% wave_atr=ANNOTD(2:end,1);
% 
% for j=1:length(splitpoint_x)
%     if j==1
%        figure(8)
%        plot(DenoisingSignal,'k');
%        axis([1,splitpoint_x(j+1),-4,4]);
%        axis off
%        print('-dpng',[location2 ,num2str(j),'.png'],'-r150');
%        
%        figure(8)
%        plot(DenoisingSignal,'k');
%        axis([splitpoint_x(j),splitpoint_x(j+2),-4,4]);
%        axis off
%        print('-dpng',[location2 ,num2str(j+1),'.png'],'-r150');
%          
%     elseif j>1 && j<length(splitpoint_x)-1
%        figure(8)
%        plot(DenoisingSignal,'k');
%        axis([splitpoint_x(j),splitpoint_x(j+2),-4,4]);
%        axis off
%        print('-dpng',[location2 ,num2str(j+1),'.png'],'-r150');
%     elseif j==length(splitpoint_x)-1
%        figure(8)
%        plot(DenoisingSignal,'k');
%        axis([splitpoint_x(j),points,-4,4]);
%        axis off
%        print('-dpng',[location2 ,num2str(j+1),'.png'],'-r150');
%     elseif j==length(splitpoint_x)
%        break;    
%     end
% end
% 
% fprintf(1,'\\n$> SAVE PICTURES FINISHED \n');





%%
%-----找出不正常波型位置-------------------------------------%
% wave_atr=ANNOTD(2:end,1);
% 
% [unnormal_wave_row,unnormal_wave_column]=find(wave_atr>1);
% 
%     for m=1:length(unnormal_wave_row)    
%     unnormal_wave_sign(m,1)=wave_atr(unnormal_wave_row(m),unnormal_wave_column(m));  
%     end
%     
%     
% [normal_wave_row,normal_wave_column]=find(wave_atr==1);
%    for n=1:length(normal_wave_row)    
%     normal_wave_sign(n,1)=wave_atr(normal_wave_row(n),normal_wave_column(n));  
%     end
    

%%
% clear x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x16 x22 x28 x31 x32 x33 x34 x37 x38
% [x2,y2]=find(unnormal_wave_sign==2);
% [x3,y3]=find(unnormal_wave_sign==3);
% [x4,y4]=find(unnormal_wave_sign==4);
% [x5,y5]=find(unnormal_wave_sign==5);
% [x6,y6]=find(unnormal_wave_sign==6);
% [x7,y7]=find(unnormal_wave_sign==7);
% [x8,y8]=find(unnormal_wave_sign==8);
% [x9,y9]=find(unnormal_wave_sign==9);
% [x10,y10]=find(unnormal_wave_sign==10);
% [x11,y11]=find(unnormal_wave_sign==11);
% [x12,y12]=find(unnormal_wave_sign==12);
% [x13,y13]=find(unnormal_wave_sign==13);
% [x14,y14]=find(unnormal_wave_sign==14);
% [x16,y16]=find(unnormal_wave_sign==16);
% [x22,y22]=find(unnormal_wave_sign==22);
% [x31,y31]=find(unnormal_wave_sign==31);
% [x32,y32]=find(unnormal_wave_sign==32);
% [x33,y33]=find(unnormal_wave_sign==33);
% [x34,y34]=find(unnormal_wave_sign==34);
% [x37,y37]=find(unnormal_wave_sign==37);
% [x38,y38]=find(unnormal_wave_sign==38);
% [x28,y28]=find(unnormal_wave_sign==28);

%%
% clear xx2 xx3 xx4 xx5 xx6 xx7 xx8 xx9 xx10 xx11 xx12 xx13 xx14 xx16 xx22 xx28 xx31 xx32 xx33 xx34 xx37 xx38
% xx2=unnormal_wave_row(x2);
% xx3=unnormal_wave_row(x3);
% xx4=unnormal_wave_row(x4);
% xx5=unnormal_wave_row(x5);
% xx6=unnormal_wave_row(x6);
% xx7=unnormal_wave_row(x7);
% xx8=unnormal_wave_row(x8);
% xx9=unnormal_wave_row(x9);
% xx10=unnormal_wave_row(x10);
% xx11=unnormal_wave_row(x11);
% xx12=unnormal_wave_row(x12);
% xx13=unnormal_wave_row(x13);
% xx14=unnormal_wave_row(x14);
% xx16=unnormal_wave_row(x16);
% xx22=unnormal_wave_row(x22);
% xx28=unnormal_wave_row(x28);
% xx31=unnormal_wave_row(x31);
% xx32=unnormal_wave_row(x32);
% xx33=unnormal_wave_row(x33);
% xx34=unnormal_wave_row(x34);
% xx37=unnormal_wave_row(x37);
% xx38=unnormal_wave_row(x38);
%%
% fid = fopen('xx8.txt','w');
% fprintf(fid,'%d OR ',xx8);
% fclose(fid);

%%

% fid = fopen('normal.txt','w');
% fprintf(fid,'%d OR ',normal_wave_row);
% fclose(fid);



%% 剪下並移動
% clc;
% clear p o
% sourcePath = 'C:\Users\AMD_107408104\Desktop\eddyhsu_paper\ECG_figure_database\ECG_10min_100(normal_unnormal)\normal';  
% targetPath = 'C:\Users\AMD_107408104\Desktop\eddyhsu_paper\TrainModel_figure\N(1)';  
% 
% for o=1:5
%     p=xx2(o,1);
%     movefile([sourcePath,'\',p,'.png'],targetPath);
% end

%% 複製並移動
% clc;
% clear p o
% sourcePath = 'C:\Users\AMD_107408104\Desktop\eddyhsu_paper\ECG_figure_database\ECG_10min_106(normal_unnormal)\unnormal\V(5)';  
% targetPath = 'C:\Users\AMD_107408104\Desktop\eddyhsu_paper\TrainModel_figure\V(5)';  
% copyfile(sourcePath,targetPath);


%% 標記所有波峰點
% [y_position,x_position,w,p]=findpeaks(DenoisingSignal);
% figure(5)
% plot(DenoisingSignal)
% string=['ECG signal ',DATAFILE,' 所有波峰標記'];
% title(string);
% grid on; axis tight;axis([2500,5500,-4,4]);
% hold on;
% plot(x_position,y_position,'o')
% legend('ECG','所有波峰標記');
% xlabel('點數');
% ylabel('電壓(mV)');

%% 標記R波
% aa=[x_position,y_position];
% [M,N]=sort(aa(:,2),'descend');
% for ii=1:length(N) 
%   bb(ii,1)=aa(N(ii,1),1);
% end
% clear N
% M=M(1:length(ANNOTD)-1,1); %降冪的R波電壓
% bb=bb(1:length(ANNOTD)-1,1);%降冪的R波電壓位置
% [M1,N1]=sort(bb,'ascend');
% for jj=1:length(bb)
%   cc(jj,1)=bb(N1(jj,1),1);
% end
% R_x=cc;
% 
% for kk=1:length(cc)
%     ccc(kk,1)=find(aa==cc(kk,1));
%     R_y(kk,1)=aa(ccc(kk,1),2);
% end
% 
% 
% figure(20)
% plot(DenoisingSignal);
% string=['ECG signal ',DATAFILE,' R點標記'];
% title(string);
% grid on; axis tight;axis([1,points,-4,4]);
% hold on;
% plot(R_x,R_y,'o');
% legend('ECG','R波');
% xlabel('點數');
% ylabel('電壓(mV)');
%% 標記P波
% x_position2=x_position;
% y_position2=y_position;
% %-----去除R波x位置-------------%
% for ll=1:length(x_position2)
%     for mm=1:length(R_x)
%        if x_position2(ll,1)==R_x(mm,1)    
%         x_position2(ll,1)=0;
%         y_position2(ll,1)=-100;
%        end
%     end
% end
% 
% for ii=1:length(ccc)
%     if ii==1
%      a1(ii,1)=max(y_position2(1:ccc(ii))); 
%      a2(ii,1)=find(y_position2==a1(ii,1));
%     elseif ii>1
%      a1(ii,1)=max(y_position2(ccc(ii)-6:ccc(ii)-1));  
%      a2(ii,1)=find(y_position2==a1(ii,1));
%     end
% end
% 
% for jj=1:length(a2)
%     a3(jj,1)=x_position2(a2(jj));
% end
% P_x=a3;
% P_y=a1;
% 
% figure(21)
% plot(DenoisingSignal);
% string=['ECG signal ',DATAFILE,' P點標記'];
% title(string);
% grid on; axis tight;axis([1,points,-4,4]);
% hold on;
% plot(P_x,P_y,'+');
% legend('ECG','P波');
% xlabel('點數');
% ylabel('電壓(mV)');
%% 標記P波起始點
% for kk=1:length(P_x)
%     a4(kk,1)=find(x_position==P_x(kk));  
% end
% for ll=1:length(a4)
%     if a4(ll)==1
%     a5(ll,1)=x_position(a4(ll));
%     elseif a4(ll)==2
%     a5(ll,1)=x_position(a4(ll));
%     elseif a4(ll)==3
%     a5(ll,1)=x_position(a4(ll));
%     else
%     a5(ll,1)=x_position(a4(ll)-3);
%     end
% end
% for mm=1:length(a5)
%     a6(mm,1)=find(x_position==a5(mm));
% end
% 
% for nn=1:length(a6)
%     P_ystart(nn,1)=y_position(a6(nn)); 
% end
% P_xstart=a5;
% 
% figure(8)
% plot(DenoisingSignal);
% string=['ECG signal ',DATAFILE,' P點起始點'];
% title(string);
% grid on; axis tight;axis([1,points,-4,4]);
% hold on;
% plot(P_xstart,P_ystart,'x');
% xlabel('點數');
% ylabel('電壓(mV)');
%% islocalmin 找局部最小值當作P波起始點
% clear mark3 flag2 Min cell
% %從P波往前找第1個局部最小值
% flag2=0;
% TF = islocalmin(DenoisingSignal);
% s_x=DenoisingSignal_x(TF);
% s_y=DenoisingSignal(TF);
% for i=1:length(P_x)
%     flag2=0;
%     for j=1:length(s_x)
%           if P_x(i,1)>s_x(j,1)
%               flag2=flag2+1;
%           elseif P_x(i,1)<s_x(j,1)
%               flag2=flag2+1;
%               mark3(i,1)=flag2;           
%               break
%           end
%     end
% end
% mark3=mark3-1;
% 
% %從第1個局部最小值往前找3個局部最小值，並找到電壓最小值的位置
% cell={};
% for l=1:length(mark3)
%     if s_x(mark3(l))<4
%     cell{l,1}(:,1)=s_x(mark3(l));
%     cell{l,1}(:,2)=s_y(mark3(l));
%     elseif s_x(mark3(l))>=4
%     cell{l,1}(:,1)=[s_x(mark3(l));s_x(mark3(l)-1);s_x(mark3(l)-2);s_x(mark3(l)-3)];
%     cell{l,1}(:,2)=[s_y(mark3(l));s_y(mark3(l)-1);s_y(mark3(l)-2);s_y(mark3(l)-3)];
%     end
% end
% for m=1:length(cell)   
%    [cell2{m,1}(:,1),cell2{m,1}(:,2)]=min(abs(cell{m,1}(:,2)));
% end
% for n=1:length(cell2)
%     Min(n,1)=cell{n,1}(cell2{n,1}(:,2),1);
%     Min(n,2)=cell{n,1}(cell2{n,1}(:,2),2);
% end
% 
% P_xstart=Min(:,1);
% P_ystart=Min(:,2);
% 
% figure(8)
% plot(DenoisingSignal);
% string=['ECG signal ',DATAFILE,' P波起始點標記'];
% title(string);
% grid on;hold on;axis tight;axis([2500,5500,-4,4]);
% plot(P_xstart,P_ystart,'x');
% legend('ECG','P波起始點');
% xlabel('點數');
% ylabel('電壓(mV)');

%%
% figure(9)
% plot(DenoisingSignal);
% string=['ECG signal ',DATAFILE,' R波標記、P波標記、P波起始點標記'];
% title(string);
% grid on; axis tight;axis([2500,5500,-4,4]);
% hold on;
% plot(R_x,R_y,'o',P_x,P_y,'+',P_xstart,P_ystart,'x');
% legend('ECG','R波','P波','P波起始點');
% xlabel('點數');
% ylabel('電壓(mV)');


%% P波起始點做分離點
%---存圖片------------------------------------%
%----圖片大小(875x656)------------------------%
%-----切割1個波形----------------------------%

% wave_atr=ANNOTD(2:end,1);
% 
% for j=12:length(wave_atr)
%     if j==1
%        figure(9)
%        plot(DenoisingSignal,'k');
%        axis([1,P_xstart(j),-4,4]);
%        axis off
%        print('-dpng',[location3 ,num2str(j),'.png'],'-r150');
%     elseif j>1&&j<length(wave_atr)
%        figure(9)
%        plot(DenoisingSignal,'k');
%        axis([P_xstart(j-1),P_xstart(j),-4,4]);
%        axis off
%        print('-dpng',[location3 ,num2str(j),'.png'],'-r150');
%     elseif j==length(wave_atr)
%        if length(P_xstart)==length(wave_atr)
%        figure(9)
%        plot(DenoisingSignal,'k');
%        axis([P_xstart(j-1),P_xstart(j),-4,4]);
%        axis off
%        print('-dpng',[location3 ,num2str(j),'.png'],'-r150');
%        
%        figure(9)
%        plot(DenoisingSignal,'k');
%        axis([P_xstart(j),points,-4,4]);
%        axis off
%        print('-dpng',[location3 ,num2str(j+1),'.png'],'-r150');
%        
%        elseif length(P_x)<length(wave_atr)
%        figure(9)
%        plot(DenoisingSignal,'k');
%        axis([P_xstart(j-1),points,-4,4]);
%        axis off
%        print('-dpng',[location3 ,num2str(j),'.png'],'-r150');   
%            
%        end
%     end
%     
% end
% 
% fprintf(1,'\\n$> SAVE PICTURES FINISHED \n');


%% 找出正常與不正常波型位置
% wave_atr=ANNOTD(2:end,1);
% %----正常-----%
% [normal_wave_row,normal_wave_column]=find(wave_atr==1);
%    for n=1:length(normal_wave_row)    
%     normal_wave_sign(n,1)=wave_atr(normal_wave_row(n),normal_wave_column(n));  
%    end
% %----不正常-----%   
% [unnormal_wave_row,unnormal_wave_column]=find(wave_atr>1);
%     for m=1:length(unnormal_wave_row)    
%     unnormal_wave_sign(m,1)=wave_atr(unnormal_wave_row(m),unnormal_wave_column(m));  
%     end
    
%%    
% figure(10)
% plot(DenoisingSignal,'k');
% axis([4085,4363,-4,4]);
% axis off






