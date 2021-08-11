%20201025
function MIT_rmd(fileName,picPath)
[ann,annType] = rdann(fileName,'atr');
ecgData = rdsamp(fileName);
avgWaveLength = round(length(ecgData)/length(ann));%平均波长
halfWL = round(avgWaveLength/2);
ecgData = Denoising(ecgData(:,1));
unnormalType = [];
% plot(ecgData());
% for i = 1:length(ann)
%     text(ann(i),ecgData(ann(i)),annType(i))
% end
for i = 1:length(ann)
    if annType(i)~= 'N'
        if annType(i) == '+'||annType(i) == '~'||annType(i) == '|'||annType(i) == '"'
            continue
        end
        if annType(i) == '/'
            annType(i) = 'm';
        end
        waveArea = ecgData(max(1,ann(i)-halfWL):min(ann(i)+halfWL,length(ecgData)),:);
        set(0,'DefaultFigureVisible', 'off')
        cd(picPath)
        [~,~] = mkdir(annType(i));
        cd(annType(i))
        [~,~] = mkdir(fileName);
        cd(fileName)
        clc;
        disp("file"+fileName+".printing:"+i+" ,type:"+annType(i))
        %     rgb = interpAndPrint(unnormalDataSplit{i},colormap(jet));
        %     rgb= imresize(rgb,[512 512]);
        plot(waveArea,'Color','black');
        axis off
        axis tight
        picName = strcat(annType(i),'_',fileName,'_',int2str(i));
        print(picName,'-djpeg','-r150');
        cd ..
        cd ..
    end
end
set(0,'DefaultFigureVisible', 'on')
end