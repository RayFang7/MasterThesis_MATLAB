
clc;clear all;close all; %#ok<CLALL>
 %#ok<*AGROW>
cd 'D:\MasterPaper\CNIS\Dataset\ptbdb\'
currentFolder = pwd;
files = dir('patient*');
names = {files.name};
patientNumber = [];
fileNumber = [];
reason = {};
AFl = {};
for i = 1 : length(names)
    cd(names{i});
    heaName = {dir('*.hea').name};
    for j = 1 : length(heaName)
        heaPath = pwd+"\"+heaName{j} ;
        fid1=fopen(heaPath,'r');
        for k = 1 : 22
            fgetl(fid1);
        end
        ROA= fgetl(fid1);
        AF= fgetl(fid1);
        fclose(fid1);
        ROA = ROA(25:length(ROA));
        AF = AF(36:length(AF));
        patientNumber = [patientNumber;names{i}];
        fileNumber = [fileNumber;heaName{j}(1:8 )];
        reason{length(reason)+1} = ROA;
        AFl{length(AFl)+1} = AF;
    end
    cd ../
end
Table2 = table(patientNumber,fileNumber,reason',AFl');
cd(currentFolder);