%20201025
function read_mark_draw(fileName,picPath)
[ann,annType] = rdann(fileName,'atr');
ecgData = rdsamp(fileName);
avgWaveLength = round(length(ecgData)/length(ann));%平均波长
halfWL = round(avgWaveLength/2);
bigHalfWL = round(halfWL*1.1);
ecgData = Denoising(ecgData());
unnormalData = zeros(0,12);
unnormalAnnIndex = [];
unnormalType = [];
unnormalDataSplit = {};
% plot(ecgData());
% for i = 1:length(ann)
%     text(ann(i),ecgData(ann(i)),annType(i))
% end
for i = 1:length(ann)
    if annType(i)~= 'N'
        waveArea = ecgData(max(1,ann(i)-halfWL):min(ann(i)+halfWL,length(ecgData)),:);
        unnormalAnnIndex = [unnormalAnnIndex;length(unnormalData)+bigHalfWL];
        unnormalData = [unnormalData;waveArea]; %#ok<*AGROW>
        unnormalDataSplit{size(unnormalDataSplit,2)+1} = waveArea; %#ok<*SAGROW,*AGROW>
        unnormalType = [unnormalType;annType(i)];
    end
end
set(0,'DefaultFigureVisible', 'off')
cd(picPath)
for i = 1:length(unnormalAnnIndex)
    [~,~] = mkdir(unnormalType(i));
    cd(unnormalType(i))
    clc;
    disp(unnormalType(i:length(unnormalAnnIndex))')
    disp(strcat ('file ',fileName,'.printing:' , int2str(i),'/',int2str(length(unnormalAnnIndex))))
    rgb = interpAndPrint(unnormalDataSplit{i},colormap(jet));
    rgb= imresize(rgb,[512 512]);
    picName = strcat(unnormalType(i),'_',fileName,'_',int2str(i));
    imwrite(rgb,picName+".jpg")
    text(unnormalAnnIndex(i),unnormalData(unnormalAnnIndex(i)),unnormalType(i))
    cd ..
end
set(0,'DefaultFigureVisible', 'on')
end