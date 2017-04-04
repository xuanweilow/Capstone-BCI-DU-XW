function [ pband_1, pband_2, pband_3, pband_4 ] = ...
           freqBandpower( data,Fs,fig )
        
       t_window = 0.02; % 20ms
       length_window = t_window/(1/Fs);
       length_end = length(data)-mod(length(data),length_window); % making sure all windows have same size
       
       i = 1;
       for n = 1:length_window:length_end
           tseries = data(n:n-1+length_window);

%            [Pxx,f] = periodogram(tseries,rectwin(length(tseries)),length(tseries),Fs); 
%            pband_tot = bandpower(Pxx,f,[0 Fs/2],'psd'); 
%            pband_1(i,1) = bandpower(Pxx,f,[1 60],'psd')/pband_tot*100;
%            pband_2(i,1) = bandpower(Pxx,f,[60 100],'psd')/pband_tot*100;
%            pband_3(i,1)= bandpower(Pxx,f,[100 300],'psd')/pband_tot*100;
%            pband_4(i,1) = bandpower(Pxx,f,[300 Fs/2],'psd')/pband_tot*100;
             
           pband_tot = bandpower(tseries,Fs,[0 Fs/2]);
           pband_1(i,1) = bandpower(tseries,Fs,[1 60])/pband_tot*100;
           pband_2(i,1) = bandpower(tseries,Fs,[60 100])/pband_tot*100;
           pband_3(i,1) = bandpower(tseries,Fs,[100 300])/pband_tot*100;
           pband_4(i,1) = bandpower(tseries,Fs,[300 Fs/2])/pband_tot*100; 
           
           i = i+1;
       end
       
       if fig == 1
           figure; plot(pband_1);hold on
           plot(pband_2);
           plot(pband_3);
           plot(pband_4);
           legend('band 1','band 2','band 3','band 4');
       end
  
end

