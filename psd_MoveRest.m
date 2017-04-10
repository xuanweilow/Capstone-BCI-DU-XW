function psd_MoveRest( data_move,data_rest,Fs )
        
       t_window = 0.02; % 20ms
       length_window = t_window/(1/Fs);
       length_end = length(data_move)-mod(length(data_move),length_window); % making sure all windows have same size
       
       figure;
       for n = 1:length_window:length_end
           series_move = data_move(n:n-1+length_window);
           series_rest = data_rest(n:n-1+length_window);
            
%            [P_move,f] = periodogram(series_move,[],length(series_move),Fs); 
%            [P_rest,f] = periodogram(series_rest,[],length(series_rest),Fs);
           
           [P_move,f] = pwelch(series_move,10,5,[1:250],Fs);
           plot(f,P_move,'r'); hold on;
           
           [P_rest,f] = pwelch(series_rest,10,5,[1:250],Fs);
           plot(f,P_rest,'b');
           
           xlabel('f');
           ylabel('PSD');
           title('PSD (20ms window)');
           legend('move','rest');
                    
       end

end