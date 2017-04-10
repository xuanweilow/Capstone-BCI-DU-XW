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
moveIndices2 = fingerMovingIndex_new(fing2,Fs,0);
% Finger 3
fing3 = fingers_nodel(:,3);
moveIndices3 = fingerMovingIndex_new(fing3,Fs,0);
% Finger 4
fing4 = fingers_nodel(:,4);
moveIndices4 = fingerMovingIndex_new(fing4,Fs,0);
% Finger 5
fing5 = fingers_nodel(:,5);
moveIndices5 = fingerMovingIndex_new(fing5,Fs,0);
% Resting Indices
restIndices = fingerRestingIndex_new(length(e1),moveIndices1,moveIndices2, ...
              moveIndices3,moveIndices4,moveIndices5);

%% PSD visualization of moving and resting 
n=10; % nth move
m=10; % mth rest
e1fing1_move1 = e1Notch_Band(moveIndices1(1,n):moveIndices1(1,n)+99);
e1_rest1 = e1Notch_Band(restIndices(1,m)+100:restIndices(1,m)+199);
psd_MoveRest(e1fing1_move1,e1_rest1,Fs);

