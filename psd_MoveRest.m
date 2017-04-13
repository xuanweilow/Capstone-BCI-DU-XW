function psd_MoveRest( data_move_before,data_move_after,data_rest,Fs )
        
       t_window = 0.02; % 20ms
       length_window = t_window/(1/Fs);
       length_end = length(data_move_after)-mod(length(data_move_after),length_window); % making sure all windows have same size
       
       figure;
       for n = 1:length_window:length_end
           series_move_before = data_move_before(n:n-1+length_window);
           series_move_after = data_move_after(n:n-1+length_window);
           series_rest = data_rest(n:n-1+length_window);
                      
           [P_move_before,f] = pwelch(series_move_before,length(series_move_before),0,[4:250],Fs);
           [P_move_after,f] = pwelch(series_move_after,length(series_move_after),0,[4:250],Fs);
           [P_rest,f] = pwelch(series_rest,length(series_rest),0,[4:250],Fs);
           
           % fft
%            L = length(series_move);
%            f = linspace(0,Fs/2,L/2+1);
%            P_move = abs(fft(series_move));
%            P_move = P_move(1:L/2+1);
%            P_rest = abs(fft(series_rest));
%            P_rest = P_rest(1:L/2+1);
           
           plot(f,P_move_before,'b'); hold on;
           plot(f,P_move_after,'k');
           plot(f,P_rest,'r');
                               
       end
       
       xlabel('f');
       ylabel('PSD');
       title('PSD (20ms window)');
       legend('move before','move after','rest');
       hold off;
end