function [output] = bpfilt(input,Fs,Fstart,Fend)

    flt = @(f)(f>=Fstart & f<=Fend); 
    
    % Symmetric Padding
    if rem(length(input),2) == 0
        in1 = input((length(input)/2):-1:2);
        in2 = input((end-1):-1:(length(input)/2));
    else
        in1 = input(((length(input)+1)/2):-1:2);
        in2 = input((end-1):-1:((length(input)+1)/2));
    end
    padd = [in1;input;in2];        
    
    % Filtering
    t = length(padd);
    ft_out = fft(padd).*flt(Fs*(0:t-1)/t)';
    output = real(ifft(ft_out));
    output = output((length(in1)+1):(length(in1)+length(input))); % Removing padded values

end