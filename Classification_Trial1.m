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
% Finger 1
fing1 = fingers_nodel(:,1);
moveIndices1 = fingerMovingIndex_new(fing1,Fs,1);
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
% Resting Indices
restIndices = fingerRestingIndex_new(length(e1),moveIndices1,moveIndices2, ...
              moveIndices3,moveIndices4,moveIndices5);
%% plot rest
t = (0:length(e1)-1)/Fs;
hold on;
restplot = zeros(length(e1),1);
for k = 1:length(restIndices(1,:))
    restplot(restIndices(1,k):restIndices(2,k))=-1;
end
plot(t,restplot,'r');  
xlabel('time (sec)');
%% PSD visualization of moving and resting 
n=10; % nth move of Finger 1
m=10; % mth rest
e1fing1_before1 = e1Notch_Band(moveIndices1(1,n)-100:moveIndices1(1,n)-1);
e1fing1_after1 = e1Notch_Band(moveIndices1(1,n):moveIndices1(1,n)+99);
e1_rest1 = e1Notch_Band(restIndices(1,m)+200:restIndices(1,m)+299);
psd_MoveRest(e1fing1_before1,e1fing1_after1,e1_rest1,Fs);
%% plot for one 20ms
n=10; % nth move of Finger 1
m=10; % mth rest
e1fing1_before1 = e1Notch_Band(moveIndices1(1,n)-20:moveIndices1(1,n)-1);
e1fing1_after1 = e1Notch_Band(moveIndices1(1,n):moveIndices1(1,n)+19);
e1_rest1 = e1Notch_Band(restIndices(1,m)+200:restIndices(1,m)+219);
psd_MoveRest(e1fing1_before1,e1fing1_after1,e1_rest1,Fs);

%% Feature extraction using Band Power
e1fing1_move1_Bands = freqBandpower(e1fing1_after1,Fs,0);    
e1_rest1_Bands = freqBandpower(e1_rest1,Fs,0);  

%% Classification using Support Vector Machine
e1fing1_move1_Band_1_3 = [e1fing1_move1_Bands(:,1), e1fing1_move1_Bands(:,3)];
e1_rest1_Band_1_3 = [e1_rest1_Bands(:,1), e1_rest1_Bands(:,3)];
% svm_train(e1fing1_move1_Band_1_3,e1fing1_rest1_Band_1_3,1,-1,1);
%% Move vs Rest Plot 
e1fing1_move1_Band_1_3_average = mean(e1fing1_move1_Band_1_3);
e1_rest1_Band_1_3_average = mean(e1_rest1_Band_1_3);

figure; 
scatter(e1fing1_move1_Band_1_3(:,1),e1fing1_move1_Band_1_3(:,2),'b'); hold on
scatter(e1_rest1_Band_1_3(:,1),e1_rest1_Band_1_3(:,2),'r');
scatter(e1fing1_move1_Band_1_3_average(1),e1fing1_move1_Band_1_3_average(2),'b','filled');
scatter(e1_rest1_Band_1_3_average(1),e1_rest1_Band_1_3_average(2),'r','filled');
xlabel('Band 1');ylabel('Band 3');
title(sprintf('%dth move of Finger 1, %dth rest',n,m));
legend('move','rest','move average','rest average');   
hold off;

%%
for n = 1:length(restIndices)
    e1Rest = e1Notch_Band(restIndices(1,n)+200:restIndices(1,n)+299);
    e1Rest_Bands = freqBandpower(e1Rest,Fs,0);
    scatter(e1Rest_Bands(:,1),e1Rest_Bands(:,3),'r'); hold on
end
for n = 1:length(moveIndices1)
    e1MoveFing1 = e1Notch_Band(moveIndices1(1,n):moveIndices1(1,n)+99);
    e1MoveFing1_Bands = freqBandpower(e1MoveFing1,Fs,0);
    scatter(e1MoveFing1_Bands(:,1),e1MoveFing1_Bands(:,3),'b');
end
%%
for n = 1:length(restIndices)
    e1Rest = e1Notch_Band(restIndices(1,n)+200:restIndices(1,n)+299);
    e1Rest_Bands = freqBandpower(e1Rest,Fs,0);
    e1Rest_Band_1_3_average(n,:) = mean([e1Rest_Bands(:,1), e1Rest_Bands(:,3)]);
end
for n = 1:length(moveIndices1)
    e1MoveFing1 = e1Notch_Band(moveIndices1(1,n):moveIndices1(1,n)+99);
    e1MoveFing1_Bands = freqBandpower(e1MoveFing1,Fs,0);
    e1MoveFing1_Band_1_3_average(n,:) = mean([e1MoveFing1_Bands(:,1), e1MoveFing1_Bands(:,3)]);
end

figure;
scatter(e1Rest_Band_1_3_average(:,1),e1Rest_Band_1_3_average(:,2)); hold on
scatter(e1MoveFing1_Band_1_3_average(:,1),e1MoveFing1_Band_1_3_average(:,2));
