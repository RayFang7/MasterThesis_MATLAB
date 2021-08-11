%clc; clear;close all;
%20201026 for ST database
%Should add read_mark_draw.m to MATLAB path before running this script.
%todo: max and min to check the R 
%single lead mode
scriptPath = pwd;
lead = 8;
filePath = "H:\Dataset\st-petersburg-incart-12-lead-arrhythmia-database-1.0.0\files\";
picPath = 'D:\MasterPaper\pics';
cd(filePath)
fileList = {dir('*.dat').name};
for fN = 1 : length(fileList)
    cd(filePath)
    read_mark_draw(erase(fileList{fN},".dat"),picPath,lead)
    disp(strcat('File ',int2str(fN),'/',int2str(length(fileList))))
end
cd(scriptPath)

%%

