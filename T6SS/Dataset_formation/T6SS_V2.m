% 相比于T6SS_V1，调整了miniBatchSize=36、Shuffle="every-epoch"、lstmLayer
clear;clc;close all;

data = readmatrix('T6SS_Positive.txt')
data2 = readmatrix('T6SS_Negative.txt')
TotalData_X = [data; data2]
a = linspace(1,1,414)
b = linspace(2,2,1111)
TotalData_Y = [a,b]
TotalData_Y = categorical(TotalData_Y)'

%% 切分训练集和测试集，70%为训练集
dataNumber = size(TotalData_X,1); %%样本个数

randIndex = randperm(dataNumber);  %%打乱数组
new_data_X = TotalData_X(randIndex,:);   
new_data_Y = TotalData_Y(randIndex,:);   

% ceil将 X 的每个元素四舍五入到大于或等于该元素的最接近整数。
testindex = ceil(dataNumber * 0.3); %% 获得分界下标,

XTest = new_data_X(1:testindex,:);
XTrain = new_data_X(testindex+1:end,:);
YTest = new_data_Y(1:testindex,:);
YTrain = new_data_Y(testindex+1:end,:);

XTrain = mat2cell(XTrain,linspace(1,1,size(TotalData_X,1)-testindex))
XTest = mat2cell(XTest,linspace(1,1,testindex))

miniBatchSize = 32

%% 定义 LSTM 网络架构
% 将输入大小指定为序列大小 1（输入数据的维度）。
% 指定具有 100 个隐含单元的双向 LSTM 层，并输出序列的最后一个元素。
% 最后，通过包含大小为2的全连接层，后跟 softmax 层和分类层，来指定2个类。

inputSize = 1;
numHiddenUnits = 100;
numClasses = 2;

layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits,OutputMode="last")
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer]

% 指定训练选项
% 指定求解器为 "adam"，梯度阈值为 1，最大轮数为 50。optimizer
% 要填充数据以使长度与最长序列相同，请将序列长度指定为 "longest"。
% 要确保数据保持按序列长度排序的状态，请指定从不打乱数据。

options = trainingOptions("adam", ...
    ExecutionEnvironment="cpu", ...
    GradientThreshold=1, ...
    MaxEpochs=50, ...
    MiniBatchSize=miniBatchSize, ...
    SequenceLength="longest", ...
    Shuffle="every-epoch", ...
    Verbose=0, ...
    Plots="training-progress");

%% 训练 LSTM 网络
% 使用trainNetwork以指定的训练选项训练 LSTM 网络。
net = trainNetwork(XTrain,YTrain,layers,options);

% 对测试数据进行分类。
% 要减少分类过程中引入的填充量，请指定使用相同的小批量大小进行训练。
% 要应用与训练数据相同的填充，请将序列长度指定为 "longest"。
YPred = classify(net,XTest, ...
    MiniBatchSize=miniBatchSize, ...
    SequenceLength="longest");

% 计算预测值的分类准确度。
acc = sum(YPred == YTest)./numel(YTest)
acc

save('T6SS_LSTM_Model.mat',net); 
