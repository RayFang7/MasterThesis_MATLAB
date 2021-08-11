clc; clear;close all;
load('Table.mat')
dbstop if error;
fileName = fileInfo.fName;
leadNo = 10;
scriptPath = pwd;
addpath("lib")
dataPath = getenv("imageDataPath");
savePath = dataPath+"\ptb12";
[~,~] = mkdir(savePath);
for fN = 1 : length(fileName) %til 144
    stringname = fileName(fN); %文件名
    cd(scriptPath)
    PTB_Load_Draw
end
set(0,'DefaultFigureVisible', 'on')