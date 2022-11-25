clear;clc;close all;

filename = 'T6SS_LSTM_Model_V4.mat';

load(filename,'-mat');

data = readmatrix('YN46_Mass.txt')

data2 = mat2cell(data,linspace(1,1,size(data,1)))

YPred = classify(net,data2, ...
    MiniBatchSize=miniBatchSize, ...
    SequenceLength="longest");