%将文件按照病人分割成三个集合
%   20201111 bulided up
%   20201112 upgraded speed
%   20201113 add split percentage
%   20201114 打乱病人顺序
%   20201115 改良折叠算法
%   20201116 真·十折
function [trainSet, devSet, testSet] = patientSplit(path,foldNo,foldTimes,percentage)
per = percentage;
oldimds = imageDatastore(path, ...
    'IncludeSubfolders',true, ...
    'LabelSource','foldernames');
files = oldimds.Files;
files = regexp( files, "patient[0-9][0-9][0-9]",'match');
files = cellfun(@(x)x(:,1), files);
pList = unique(files);
pLen = length(pList);
rng('default')
pList = pList(randperm(pLen)); 
rng('shuffle')
for i = 0:foldTimes-1
    up = round((1/foldTimes)*i*pLen)+1;
    down = round((1/foldTimes)*pLen*(i+1));
    index = [];
    for j = up : down
        index = [index;find(contains(files,pList(j)))]; %#ok<AGROW>
    end
    subds{i+1} = index; %#ok<AGROW>
end
indexList = [1:foldTimes 1:foldTimes];
% testIndex = foldNo;
% devIndex = [indexList(foldNo+1) indexList(foldNo+2)];
% trainIndex = indexList(foldNo+3:foldNo+foldTimes-1);
% trainSet = subdsCombine(subds(trainIndex), oldimds , per);
% devSet = subdsCombine(subds(devIndex), oldimds , per);
% testSet = subdsCombine(subds(testIndex), oldimds , 1);
%% real ten fold
devIndex = indexList(foldNo);
trainIndex = indexList(foldNo+1:foldNo+foldTimes-1);
trainSet = subdsCombine(subds(trainIndex), oldimds , per);
devSet = subdsCombine(subds(devIndex), oldimds , per);
testSet = {};
end


function imdOutput = subdsCombine(cellInput,imdInput , per)
index = cellInput{1};
index = index(randperm(length(index)));
index = index(1:round(length(index)/(1/per)));
for i = 2:length(cellInput)
    temp = cellInput{i};
    temp = temp(randperm(length(temp)));
    temp = temp(1:round(length(temp)/(1/per)));
    index = [index;temp]; %#ok<AGROW>
end
imdOutput = subset(imdInput,index);
end
