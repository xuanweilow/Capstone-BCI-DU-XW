function output = testBCI(data,Fs,Spfilt,lda_W,lda_B) 

    % Spatial filtered
    data_filt = data*Spfilt;
    % Bandpass filtered
    for i = 1:size(data_filt,2)
        data_filt(:,i) = bpfilt(data_filt(:,i),Fs,7,30);
    end

    % Classifier output
    output = log(var(data_filt))*lda_W - lda_B;
    if sign(output)<0
        output = 1;
    else
        output = 2;
    end
    
end