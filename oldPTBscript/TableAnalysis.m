clc;clear all;
load('Table.mat')
row = (Reasonofadmission.Reason == "Healthy control");
healthyName = Reasonofadmission.fileName(row);
row = (Reasonofadmission.fileName == healthyName(5));
Reasonofadmission.patientName(row)  