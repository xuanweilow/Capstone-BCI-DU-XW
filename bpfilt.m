function [output] = bpfilt(input,Fs,Fstart,Fend)
    % input_c takes a matrix whose rows are a time series data & columns are variables

    flt = @(f)(f>=Fstart & f<=Fend); 
%     flt = designfilt('bandpassfir','StopbandFrequency1',Fstart-1,'PassbandFrequency1',Fstart,'PassbandFrequency2',Fend,...
%                      'StopbandFrequency2',Fend+1,'StopbandAttenuation1',50,'PassbandRipple',1,'StopbandAttenuation2',50,...
%                      'DesignMethod','kaiserwin','SampleRate',Fs);
% %     fvtool(flt); for filter analyse
    
    for c = 1:size(input,2)
        input_c = input(:,c);
        
        % Symmetric Padding
        if rem(length(input_c),2) == 0
            in1 = input_c((length(input_c)/2):-1:2);
            in2 = input_c((end-1):-1:(length(input_c)/2));
        else
            in1 = input_c(((length(input_c)+1)/2):-1:2);
            in2 = input_c((end-1):-1:((length(input_c)+1)/2));
        end
        padd = [in1;input_c;in2];        

        % Filtering
        t = length(padd);
        ft_out = fft(padd).*flt(Fs*(0:t-1)/t)';
        outpad_c = real(ifft(ft_out));
%         outpad_c = filtfilt(flt,padd);

        output(:,c) = outpad_c((length(in1)+1):(length(in1)+length(input_c))); % Removing padded values
    end

end