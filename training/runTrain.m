%运行训练
%   20201114 增加随机处理
function [net,options,traininfo,result] = runTrain(params,DatasetPath)
[layers,params.picSize] = getLayers(params.netName);
[trainSet, devSet, testSet, numClasses,classWeights] = getDataset(params,DatasetPath); %#ok<ASGLU>
tempLayers = [
    fullyConnectedLayer(numClasses,"Name","fullyConnectedLayer")
    softmaxLayer("Name","softmaxLayer")
    classificationLayer("Name","ClassificationLayer_predictions")];
%     labelSmoothClassificationLayer("SmoothClassoutput",numClasses,params.labelSmoothFactor)];
%     weightedClassificationLayer("WeightedClassoutput",classWeights)];
layers = addLayers(layers,tempLayers);
layers = connectLayers(layers,"output","fullyConnectedLayer");
clear tempLayers
% layers = resnetV218(numClasses,classWeights);
disp("using "+ params.netName + " network")
options = getOptions(params,devSet);

[net, traininfo] = trainNetwork(trainSet,layers,options);
result = Mclassify(net,params.picSize,1,testSet);
meanResult = resultCalc({result},{traininfo});
disp("preSum = " + meanResult.preSum...
    +" reSum = " + meanResult.reSum...
    +" f1Sum = " + meanResult.f1Sum...
    +" valSum = " + meanResult.valSum)
end
%% getLayers
function [layers,picSize] = getLayers(netName)
switch (netName)
    case 'alexNet'
        layers = MalexNet();
        picSize = [227 227];
    case 'googlenet'
        layers = MgoogLenet();
        picSize = [224 224];
    case 'resnet18'
        layers = Mresnet18();
        picSize = [224 224];
    case 'resnet50'
        layers = resnet50('Weights','none');
        layers = removeLayers(layers,'fc1000'); %resnet
        layers = removeLayers(layers,'ClassificationLayer_fc1000'); %resnet18
        picSize = [224 224];
    case 'resnetSmall'
        layers = resnetSmall();
        picSize = [224 224];
end
end

%% getOptions
function options = getOptions(params,devSet)
options = trainingOptions(params.gdMethod, ...
    'MaxEpochs',params.MaxEpochs, ...
    'ValidationData',devSet, ...
    'ValidationFrequency',15, ...
    'InitialLearnRate',params.InitialLearnRate,...
    "MiniBatchSize",params.MiniBatchSize,...
    'ValidationPatience',20, "Plots","training-progress",...
    "L2Regularization",params.L2Re,...
    "LearnRateDropFactor",0.45,...
    "LearnRateDropPeriod",1,...
    "LearnRateSchedule","piecewise",...
    "Shuffle","every-epoch");
disp('Options created');
end

%% getDataset
function [trainSet, devSet, testSet,numClasses,classWeights] = getDataset(params,DatasetPath)
[trainSet, devSet, testSet] = patientSplit(DatasetPath+"\H",params.foldNo,params.foldTimes,params.dataPercent);
[trainSet2, devSet2, testSet2] = patientSplit(DatasetPath+"\MI",params.foldNo,params.foldTimes,params.dataPercent);
trainSet = imdCat(trainSet,trainSet2);
devSet = imdCat(devSet,devSet2);
% testSet = imdCat(testSet,testSet2);

trainSet = shuffle(trainSet);
devSet = shuffle(devSet);
% testSet = shuffle(testSet);
testSet  = devSet; %%for result calc
imageAugmenter = imageDataAugmenter( ...
    'RandXTranslation',[-20 20], ...
    'RandYTranslation',[-20 20], ...
    'RandXScale', [0.7,1.3],...
    'RandYScale', [0.7,1.3]);

categ = categories(trainSet.Labels);
fprintf("-------------------------------------------------------------------\n");
disp   ("                              Fold no."+params.foldNo);
fprintf("%20s             | %9s| %9s| %6s\n","Label","Amount","Percentage","Weight")
clear oldimds;
numClasses = numel(categories(trainSet.Labels));
for i = 1:numClasses
    classWeights(i) = length(trainSet.Files)/sum(trainSet.Labels==categ{i}); %#ok<AGROW>
end
classWeights = classWeights/(max(classWeights));
for i = 1:length(categ)
    fprintf("   %-30s|%10d|%10.2f%%|%8.4f \n",categ{i},sum(trainSet.Labels==categ{i}),sum(trainSet.Labels==categ{i})/length(trainSet.Files)*100,classWeights(i));
end
trainSet = augmentedImageDatastore(params.picSize,trainSet,'DataAugmentation',imageAugmenter);
devSet = augmentedImageDatastore(params.picSize,devSet);
disp('Data preprocess finished');

%silly but steady code to make sure there isn't data spoil.
setCombine(1:length(devSet.Files)) = devSet.Files;
% setCombine(end+1:end+length(testSet.Files)) = testSet.Files;
setCombine(end+1:end+length(trainSet.Files)) = trainSet.Files;
devPat = getPatient(devSet.Files);
trainPat = getPatient(trainSet.Files);
% testPat =  getPatient(testSet.Files);

patCombine(1:length(devPat)) = devPat;
patCombine(end+1:end+length(trainPat)) = trainPat;
% patCombine(end+1:end+length(testPat)) = testPat;

if length(unique(setCombine))==length(setCombine) && length(unique(patCombine))==length(patCombine)
    disp('Data safty checked')
end

end

%% functions
function output = getPatient(cellin)
files = regexp( cellin, "patient[0-9][0-9][0-9]",'match');
files = cellfun(@(x)x(:,1), files);
output = unique(files);
end

function imd = imdCat(imd1,imd2)
imdTemp = imageDatastore(cat(1,imd1.Files,imd2.Files));
imdTemp.Labels = (cat(1,imd1.Labels,imd2.Labels));
imd = imdTemp;
end
