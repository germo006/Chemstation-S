%% Noah Germolus 18 Oct 2020
% This is more or less a driver for ConcatData.m

clear; clc
listFiles = ['AE2123_NPG1.2022.03.01.mat';...
             'AE2123_NPG2.2022.03.01.mat'];
listFiles=string(listFiles);

ConcatData_v1(listFiles, 'AE2123_BC_ZippedC12.mat')
clear listFiles