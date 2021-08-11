%20201031 Automatic label detect version
%result : prob
%#ok<*AGROW>
currentFolder = pwd;
path = "D:\MasterPaper\CNIS\Workspace\Data\202\test2";
network = trainedNetwork_5;
cd(path);
folders = {dir('*').name};
prob = cell(2,length(folders)-2);
for i = 3:length(folders)
    prob{1,i-2} = folders{i};
    cd(folders{i})
    files = dir('*');
    names = {files.name};
    result = cell(length(names)-2,2);
    for j = 3 : length(names)
        I = imread(names{j});
        I = imresize(I, [224 224]);
        [YPred,probs] = classify(network,I);
        result{j-2,1}=probs;
        result{j-2,2}=string(YPred);
    end
    prob{2,i-2}=result;
    cd ..
end
cd(currentFolder);