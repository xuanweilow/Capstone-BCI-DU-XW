%% Main page
%% Training:
% mrk_train tags the class label at the start of each trial 
Fs = 100;
wnd_test = 0.5*Fs; % classification length (sec *Fs)
[Spfilt,lda_W,lda_B0,means_k,train_epocAdj] = trainBCI(eTrain,lblTrain,wnd_test,Fs,1);

%% Online Testing:
out = [];

% OPTION:
paramU = 0;
% assumes that 'eTest' is in chunk of 'wnd' length
features_tst{1} = []; features_tst{2} = []; % store features for later analysis
if lblTest ~= 0 % current class label
    [out(i),fea,means_k] = testBCI_1(eTest,Fs,Spfilt,lda_W,means_k,paramU);
    features_tst{lblTest} = vertcat(features_tst{lblTest},fea);
else 
    [out(i),~,~] = testBCI_1(eTest,Fs,Spfilt,lda_W,means_k,paramU);
end

%% (Analysis: test classification output visualise)
pca_data = pca([features_tst{1};features_tst{2}]); 
for k = 1:2
    x_axs{k} = features_tst{k}*lda_W;
    y_axs{k} = features_tst{k}*pca_data(:,1);
end
figure; hold on
scatter(x_axs{1},y_axs{1},'r','filled');
scatter(x_axs{2},y_axs{2},'b','filled'); 
yy = ylim; line([lda_B0 lda_B0],[yy(1) yy(2)]); 

%%
clear;clc;%close all;   
load data_set_IVa_aw.mat
load true_labels_aw.mat
Fs = 100;
E = [34 38 51 53 55 57 70 74]; % Chosen electrodes
eTrain = 0.1*double(cnt(1:mrk.pos(57)-1,E)); % train data
wnd_test = 0.5*Fs; %%
mrk_train = horzcat(mrk.pos(1:56)',true_y(1:56)');
trial_len = 3.5*Fs;
train_epoc = [0:trial_len-1];
lblTrain = zeros(size(eTrain,1),1); % labelling (0,1,2)
for j = 1:size(mrk_train,1)
   lblTrain(mrk_train(j,1)+train_epoc) = mrk_train(j,2);
end

