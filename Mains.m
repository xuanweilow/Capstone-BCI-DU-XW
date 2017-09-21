%% Main page
%% Training:
Fs = 250;
wnd = 0.5; % classification length (sec)
[train_data,avrC,Spfilt,lda_W,lda_B] = trainBCI(train_data,mrk,wnd,Fs,1); 

%% Online Testing:
out = [];
% assumes that 'test_data' is in chunk of 'wnd' length
out(i) = testBCI(test_data,Fs,Spfilt,lda_W,lda_B);