clear;close all;clc;
% BCL Competition Data Set 4/sub1 
load('sub1_comp.mat');
% %Plotting Finger Movements
% fingers = cellstr(['thumb '; 'index '; 'middle'; 'ring  '; 'little']); %Plot of 37ms and normalized Finger Positions
% figure
% for n = 1:5
%     suBandlot(5,1,n);plot(train_dg(:,n));ylabel(fingers(n));
% end
%% Filtering 
% Filtering the signal from the first electrode with a Chebyshev2 Bandpass 
% filter with a 4-250Hz pass-band and a 60Hz Notching filter
Fs = 1000; tn = linspace(0,399.999,400000)';
e1 = train_data(:,1);
% [e1Band,e1FBand] = FilteringT(e1,tn,Fs,'cheb2',01);%Bandpass Filtering with plots
% [e1Band_Notch, e1FBand_Notch] = FilteringT(e1Band,tn,Fs,'notch',01);%Notch Filtering with plots
[e1Notch,e1FNotch] = FilteringT(e1,tn,Fs,'notch',00);%Notch Filtering with plots
[e1Notch_Band, e1FNotch_Band] = FilteringT(e1Notch,tn,Fs,'cheb2',00);%Bandpass Filtering with plots
%% Removing Finger Delay (37ms)
fingers_nodel = [train_dg(37:end,:); train_dg((end-35):end,:)];
%% Obtaining the moving and resting regions of the hand
fing1 = fingers_nodel(:,1);
moveIndices1 = fingerMovingIndex_new(fing1,Fs,1);
restIndices1 = fingerRestingIndex_new(fing1,moveIndices1);
% for n = 1:length(movingIndices1(1,:))
%     e1_move_n = e1Notch_Band(moveIndices1(1,n):moveIndices1(2,n)); ...
%     
% end
%% 
% Finger 2
fing2 = fingers_nodel(:,2);
moveIndices2 = fingerMovingIndex_new(fing2,Fs,1);
% Finger 3
fing3 = fingers_nodel(:,3);
moveIndices3 = fingerMovingIndex_new(fing3,Fs,1);
% Finger 4
fing4 = fingers_nodel(:,4);
moveIndices4 = fingerMovingIndex_new(fing4,Fs,1);
% Finger 5
fing5 = fingers_nodel(:,5);
moveIndices5 = fingerMovingIndex_new(fing5,Fs,1);

%% Feature extraction using Band Power
n=5; % nth move
e1fing1_move1 = e1Notch_Band(moveIndices1(1,n):moveIndices1(2,n)); 
e1fing1_rest1 = e1Notch_Band(restIndices1(1,n):restIndices1(2,n));
% e1fing1_move1 = e1Notch_Band(moveIndices1(1,n):moveIndices1(1,n)+199); 
% e1fing1_rest1 = e1Notch_Band(restIndices1(2,n)-199:restIndices1(2,n));

[e1fing1_move1_Band1,e1fing1_move1_Band2,e1fing1_move1_Band3,e1fing1_move1_Band4] = freqBandpower(e1fing1_move1,Fs,1);    
[e1fing1_rest1_Band1,e1fing1_rest1_Band2,e1fing1_rest1_Band3,e1fing1_rest1_Band4] = freqBandpower(e1fing1_rest1,Fs,1);  

e1fing1_move1_Band1_average = mean(e1fing1_move1_Band1);
e1fing1_move1_Band3_average = mean(e1fing1_move1_Band3);
e1fing1_rest1_Band1_average = mean(e1fing1_rest1_Band1);
e1fing1_rest1_Band3_average = mean(e1fing1_rest1_Band3);

%% Classification using Support Vector Machine
e1fing1_move1_Bands = [e1fing1_move1_Band1, e1fing1_move1_Band3];
e1fing1_rest1_Bands = [e1fing1_rest1_Band1, e1fing1_rest1_Band3];
svm_train(e1fing1_move1_Bands,e1fing1_rest1_Bands,1,-1,1);
%%
figure; 
scatter(e1fing1_move1_Bands(:,1),e1fing1_move1_Bands(:,2),'b'); hold on
scatter(e1fing1_rest1_Bands(:,1),e1fing1_rest1_Bands(:,2),'r');
xlabel('Band 1');ylabel('Band 3');
plotaxis = 0:100;
plot(e1fing1_move1_Band1_average.*ones(size(plotaxis)),plotaxis);
plot(plotaxis,e1fing1_move1_Band3_average.*ones(size(plotaxis)));
plot(e1fing1_rest1_Band1_average.*ones(size(plotaxis)),plotaxis);
plot(plotaxis,e1fing1_rest1_Band3_average.*ones(size(plotaxis)));
legend('move','rest','move band1 average','move band3 average','rest band1 average','rest band3 average');   
