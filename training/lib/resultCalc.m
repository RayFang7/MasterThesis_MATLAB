function [meanResult] = resultCalc(result,traininfo)
preSum = 0;
reSum = 0;
f1Sum = 0;
valSum = 0;
for i = 1:length(result)
    preSum = preSum+result{i}.precision;
    reSum = reSum + result{i}.recall;
    f1Sum = f1Sum + result{i}.meanF1;
    valSum = valSum + traininfo{i}.FinalValidationAccuracy;
end
meanResult.preSum = mean(preSum)/length(result);
meanResult.reSum = mean(reSum)/length(result);
meanResult.f1Sum = mean(f1Sum)/length(result);
meanResult.valSum = mean(valSum)/length(traininfo)/100;
end