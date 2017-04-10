function [ restIndice ] = fingerRestingIndex_new( sample_length,fing1_moveIndice, ...
                                                           fing2_moveIndice, ...  
                                                           fing3_moveIndice, ...
                                                           fing4_moveIndice, ...
                                                           fing5_moveIndice )

    fing1_move = zeros(sample_length,1);
    for k = 1:length(fing1_moveIndice(1,:))
        fing1_move(fing1_moveIndice(1,k):fing1_moveIndice(2,k))=1;
    end
    fing2_move = zeros(sample_length,1);
    for k = 1:length(fing2_moveIndice(1,:))
        fing2_move(fing2_moveIndice(1,k):fing2_moveIndice(2,k))=1;
    end
    fing3_move = zeros(sample_length,1);
    for k = 1:length(fing3_moveIndice(1,:))
        fing3_move(fing3_moveIndice(1,k):fing3_moveIndice(2,k))=1;
    end
    fing4_move = zeros(sample_length,1);
    for k = 1:length(fing4_moveIndice(1,:))
        fing4_move(fing4_moveIndice(1,k):fing4_moveIndice(2,k))=1;
    end
    fing5_move = zeros(sample_length,1);
    for k = 1:length(fing5_moveIndice(1,:))
        fing5_move(fing5_moveIndice(1,k):fing5_moveIndice(2,k))=1;
    end
    
    min_rest_length = 1000;
    
    new_rest = 1;
    i = 1;
    for n = 1:sample_length
        if fing1_move(n) + fing2_move(n) + fing3_move(n) + fing4_move(n)...
           + fing5_move(n) == 0
            if new_rest == 1
                rest_start = n;
                new_rest = 0;
            end
        else
            if new_rest == 0 
                if n-rest_start >= min_rest_length
                    rest_end = n-1;
                    restIndice(:,i) = [rest_start rest_end];
                    i = i+1;
                end
                new_rest = 1;
            end
        end
    end        
    
end

