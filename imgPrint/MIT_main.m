%clc; clear;close all;
%20201107 for MITDB
scriptPath = pwd;
filePath = "D:\MasterPaper\eddy_paper_data\MIT_BIH_database";
picPath = 'D:\MasterPaper\MITpics';
cd(filePath)
fileList = {dir('*.dat').name};
for fN = 35 : length(fileList)
    cd(filePath)
    fileName = erase(fileList{fN},".dat");
    MIT_rmd(fileName,picPath)
    disp(strcat('File ',int2str(fN),'/',int2str(length(fileList))))
end
cd(scriptPath)z

%%

