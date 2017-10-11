%% Main page
%% Training:
% mrk_train tags the class label at the start of each trial 
Fs = 250;
wnd_test = 0.5*Fs; % classification length (sec *Fs)
[Spfilt,lda_W,lda_B0,means_k,train_epocAdj] = trainBCI(eTrain,lblTrain,wnd_test,Fs,1);

%% Online Testing:
% out = []; 
% lblTest is the current class label
% OPTION:
paramU = 0;
% assumes that 'eTest' is in chunk of 'wnd' length
features_tst{1} = []; features_tst{2} = []; % store features for later analysis
if lblTest ~= 0 % current class label
    [out,fea,means_k] = testBCI(eTest,Fs,Spfilt,lda_W,means_k,paramU);
    features_tst{lblTest} = vertcat(features_tst{lblTest},fea);
else 
    [out,~,~] = testBCI_1(eTest,Fs,Spfilt,lda_W,means_k,paramU);
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
