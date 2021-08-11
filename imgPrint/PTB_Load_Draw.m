%20201030 I hope to add version control system.
%20201102 Changed print method to interp and print
%20201104 Updated R finding method
%20201108 AF
%20201114 改为全导程
%%
global sfreq
set(0,'DefaultFigureVisible', 'off')
load('cm.mat');
cm = colormap(jet);
str = fN+"/516,"+round(fN/516*100,2)+"%";
row = (fileInfo.fName == stringname);
patientNumber = fileInfo.pNumber(row)  ;
ROA = fileInfo.ROA(row);
AF = fileInfo.AF(row);
PATH= "D:\MasterPaper\CNIS\Dataset\ptbdb\"+patientNumber+"\"; % 储存路径
HEADERFILE= strcat(stringname,'.hea');     % header-file in text format
ATRFILE= strcat(stringname,'.atr');        % attributes-file in binary format
DATAFILE=strcat(stringname,'.dat');        % data-file
XYZFILE=strcat(stringname,'.xyz');        % data-file
leadname=["I","II","III","avr","avl","avf","v1","v2","v3","v4","v5","v6","vx","vy","yz"];
%123456/L3 AVF L2 -avr L1 AVL
%------ LOAD HEADER DATA --------------------------------------------------
signalh= fullfile(PATH, HEADERFILE);
fid1=fopen(signalh,'r');
z= fgetl(fid1);
A= sscanf(z, '%*s %d %d %d',[1,3]);
points = A(3);
sfreq=A(2);   % sample rate of data取样频率
clear A;
for k=1:12
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
M = loadfile(fullfile(PATH, DATAFILE),12,points);
for i = 1:12
    M(:,i)=(M(:,i)-zerovalue(i))/gain(i);
end
M=[M(:,3) M(:,6) M(:,2) -M(:,4) M(:,1) M(:,5) M(:,7:12)] ; %重新排列;
%M=M(2000:25000,:);%提取测试资料
points = length(M);
if ROA == "MI"
    M = M(1:min(38400,points),:);
end
TIME = (0:(points-1))/sfreq;
%-------------------------------------------------
M=Denoising2(M);


% Rerr=0;
% try
%     F_FindR
% catch
%     disp(fN+" error")
%     clearvars -except M te points sfreq stringname leadname Rloc picName fileName fN Reasonofadmission wb rw Rerr
%     return
% end
Rloc = findR(M);  %cite :Complete Pan Tompkins Implementation ECG QRS dete
%%
for a = 1:length(Rloc)
    plot(Rloc(a),M(Rloc(a)),'o','color','g')
end
CutLoc = Rloc-199;
CutEnd = Rloc+400;
%%
%绘图
%figure(1)
%[M1,c]=contour(1:points,1:15,M(1:points,:)',100); %显示所有值
%R值取值确认
%line(repelem(Rloc,15,1),1:15,"Color","white","LineWidth",1,'LineStyle','- -');
%line(repelem(CutLoc,15,1),1:15,"Color","blue","LineWidth",1,'LineStyle','- -');
% pList = ["hp";"fp";"ap"];
for i = 3:length(CutLoc)-1
    cd(savePath)
    pName = "ap";
    subPath = pName+"/"+ROA+"/"+AF;
    [~,~] = mkdir(subPath);
    cd(subPath)
    left = max(1,round(CutLoc(i)));
    right = min(points,round(CutEnd(i)));
    I2 = interpAndPrint(M(left:right,:),cm);
    I2 = imresize(I2,[512 512]);
    picName=ROA+"_"+pName+"_"+AF+"-"+patientNumber+"_"+stringname;
    imwrite(I2,picName+i+'.jpg');
    clc
    fprintf('file no.%d/%d \nwave no.%d/%d',fN,length(fileName),i,length(CutLoc)-1)
end
% cd ..
cd ../
%%
function data=loadfile(filepath,nosig,points)
fid2=fopen(filepath,'r');
data= fread(fid2, [nosig, points], 'uint16')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);                                 % 矩阵有3行，每行8位，= 2 * 12位
for i=1:nosig
    PR = bitand(data(:,i),32768);
    data(:,i)=bitand(data(:,i),32767)-PR;
end
end