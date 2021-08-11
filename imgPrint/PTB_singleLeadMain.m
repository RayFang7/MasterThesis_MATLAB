%20201025
%20201104 Upgraded R peak search algorithum
%Should add read_mark_draw.m to MATLAB path before running this script.
%single lead mode
clc; clear;close all;
scriptPath = pwd;
load('Table.mat')
fileName = fileInfo.fName;
addpath('lib')
picPath = 'D:\MasterPaper\PTBpics';
for i = 5 : 12
for fN = 371 : length(fileName)
    row = (fileInfo.fName == fileName(fN));
    patientNumber = fileInfo.pNumber(row)  ;
    ROA = fileInfo.ROA(row);
    filePath = "D:\MasterPaper\CNIS\Dataset\ptbdb\"+patientNumber+"\";
    cd(filePath)
    try
    PTB_rmd(fileName(fN),picPath,i,ROA,patientNumber)
    catch
    end
    disp(strcat('File ',int2str(fN),'/',int2str(length(fileName(fN)))))
end
end
cd(scriptPath)
