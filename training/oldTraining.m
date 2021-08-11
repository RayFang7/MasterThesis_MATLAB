clc;clear all
%%
params.WeightLearnRateFactor = 2;
params.BiasLearnRateFactor = 2;
params.myInitialLearnRate = 0.001;
params.optionType = 3;
params.foldNo = 1;
params.L2Re = 0.01;
i = 1;
digitDatasetPath = "D:\MasterPaper\CNIS\Workspace\PTB\pics\ap\MI";
%%
for i = 1:3
    params.foldNo = i;
    [net{i},options{i},traininfo{i},result{i}]  = uniTrain(params,digitDatasetPath);
end
%%
digitDatasetPath = "D:\MasterPaper\TrainingData\306\hp\TRAIN";
for i = 11:13
    params.foldNo = i-10;
    [net{i},options{i},traininfo{i},result{i}]  = uniTrain(params,digitDatasetPath);
end
%%
sum =0;
sum2 =0;
for i = 11:13
    sum=sum+mean(result{i}.precision);
    sum2=sum2+length(traininfo{i}.TrainingLoss);
end