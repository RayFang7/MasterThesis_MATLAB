I = imread("C:\MasterPaper\Dataset\DataVersion4\AB\B\Myocardial infarction_patient004_s0020bre_3.jpg");
I = imresize(I, [224 224]);
[YPred,probs] = classify(trainedNetwork_1,I);
imshow(I)
label = YPred;
title(string(label) + ", " + num2str(100*max(probs),3) + "%");