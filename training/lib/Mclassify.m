%classify
%I am an idiot
%type 1 for ds input, other for path input
%todo 补完分类
function result = Mclassify(network,picSize,type,Input)
if type == 1
    imds = Input;
else
imds = imageDatastore(Input, ...
    'IncludeSubfolders',true, ...
    'LabelSource','foldernames');
end
testImds = augmentedImageDatastore(picSize,imds, ...
    'ColorPreprocessing','gray2rgb');
[YPred,~] = classify(network,testImds);
if (imds.Labels(1)== "H"||imds.Labels(1)== "MI")
confusionMat = confusionmat(imds.Labels,YPred,'order',{'MI','H'});
else
confusionMat = confusionmat(imds.Labels,YPred);
end
result.confusionMat = confusionMat;
result.precision = diag(confusionMat)./sum(confusionMat,2);
result.recall = diag(confusionMat)./sum(confusionMat,1)';
result.f1Scores = 2*(result.precision.*result.recall)./(result.precision+result.recall);
result.meanF1 = mean(result.f1Scores);
end