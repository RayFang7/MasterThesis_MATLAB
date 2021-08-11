%20201025 Single lead ptb rmd
%20201104 upgrade Rpeak
function PTB_rmd(fileName,picPath,lead,ROA,patientNumber)
fileName = convertStringsToChars(fileName);
ecgData = rdsamp(fileName);
try
ecgData = Denoising(ecgData(:,13:15));
catch
    return
end


if ROA == "MI"
    ecgData = ecgData(1:20000);
end
Rmax = findR(ecgData);
avgWaveLength = round(length(ecgData)/length(Rmax));
halfWL = round(avgWaveLength/2);
bigHalfWL = round(halfWL*1.1);
set(0,'DefaultFigureVisible', 'off')
[~,~] = mkdir(picPath+"\"+lead+"\"+ROA);
cd(picPath+"\"+lead+"\"+ROA)
for i = 1:length(Rmax)
    waveArea = ecgData(max(1,Rmax(i)-bigHalfWL):min(length(ecgData),Rmax(i)+bigHalfWL),:);
    clc;
    waveArea = normalize(waveArea);
    disp(strcat ('file ',fileName,'.printing:' , int2str(i),'/',int2str(length(Rmax))))
    plot(waveArea,'Color','black');
    axis off
    axis tight
    print(strcat(ROA,'_L',int2str(lead),'_',patientNumber,'_',fileName,'_',int2str(i)),'-djpeg','-r150')
end
cd ..
set(0,'DefaultFigureVisible', 'on')
end