% 新版训练模块
%   20201111 组合了训练脚本
%   20201113 拆分Layer的生成
%   20201114 再次拆分
%#ok<*SAGROW>
%% 初始化
reset(gpuDevice(1))
clc;clearvars -except result params options traininfo netLength;
scriptPath = pwd;
cd (scriptPath)
addpath("lib")
addpath("layers")
addpath("netWork")
dataPath = getenv("imageDataPath");
%% 
params.experimentNo = 404;
params.InitialLearnRate = 0.001;
params.foldNo = 1;
params.foldTimes = 20;
params.L2Re = 0.05;
params.gdMethod = 'sgdm';
params.labelSmoothFactor = 0;
params.netName = 'alexNet';
params.dataPercent = 0.3;
params.MaxEpochs = 10;
params.MiniBatchSize = 128;
DatasetPath = dataPath + "\308\ap\TRAIN";
% DatasetPath = "D:\MasterPaper\PTBpics\4\TRAIN"
%%
if exist('netW','var') == 1 
    netLength = length(netW);
else
    netLength = 0;
end
for i = 1:20
    params.foldNo = i;
    DatasetPath = dataPath + "\308\ap\TRAIN";
    [netW{netLength+1},options{netLength+1},traininfo{netLength+1},result{netLength+1}] = runTrain(params,DatasetPath);
    paramsBackup{netLength+1} = params; 
    netLength = length(netW);
end

%%
% 
sur.val = 0;
sur.cm = [0 0;0 0];
sur.pre = [0 ;0];
sur.rec = [0;0];
sur.f1 = [0;0];
for i = 1:length(result)
    sur.val(i) = traininfo{i}.ValidationAccuracy(end);
    sur.cm = sur.cm + result{i}.confusionMat;
    sur.pre =  sur.pre  +  result{i}.precision/length(result);
    sur.rec = sur.rec + result{i}.recall/length(result);
    sur.f1 = sur.f1 +  result{i}.meanF1/length(result);
end

%%
a=devSet.files
for i = 1:length(A)
name = A{i};
row = (fileInfo.pNumber == name);
AFc{i} = unique(fileInfo.AF(row));
end