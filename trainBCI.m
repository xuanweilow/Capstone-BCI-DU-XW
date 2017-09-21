function [train_data,avrC,Spfilt,lda_W,lda_B] = trainBCI( eTrain,mrk,wnd,Fs,seeplot ) 
% (note: mrk in sparse, labels 1&2, assuming class data equal and div by 10)

    % OPTIONS:
    freqRange = [7 30]; 
    wnd_train = [1 3]; 
    nof = 4; % needs validation  

    %% Bandpass filtering
    eTrain_flt = zeros(size(eTrain));
    for i = 1:size(eTrain,2)
        eTrain_flt(:,i) = bpfilt(eTrain(:,i),Fs,freqRange(1),freqRange(2));
    end
    eTrain = eTrain_flt; 

    %% Epoching 
%     %(assumes rest class between trials)
%     wnd_test = wnd*Fs;
%     nstep = wnd_test/2;
%     train_epoc = [round(Fs*wnd_train(1)) round(Fs*wnd_train(2))];     
%     epoc_idx = train_epoc(1):nstep:train_epoc(2); % make it divisible for wnd_test purpose
%     if mod(train_epoc(2)-train_epoc(1)+1,nstep) ~= 0
%        epoc_idx = epoc_idx(1:end-1);
%     end
%     nclass{1} = 0; nclass{2} = 0; % number of data trials from class 
%     pos = 1;
%     for i = 2:length(label)
%         pos = pos+1;
%         k = label(i);
%         if k~= 0 && k~= label(i-1) 
%             pos = 1;
%         end
%         if k~= 0 && nnz(pos == epoc_idx)~= 0 
%             nclass{k} = nclass{k}+1;
%             tmp_data = eTrain(pos:pos+wnd_test-1,:);
%             train_data{k,nclass{k}} = tmp_data; % store train epochs
%             tmp_C = tmp_data'*tmp_data; % sample covariance
%             C{k}(:,:,nclass{k}) = tmp_C./trace(tmp_C);
%         end
%     end
%     avrC{1} = mean(C{1},3); % class average covariance
%     avrC{2} = mean(C{2},3);
    
    wnd_test = wnd*Fs;
    train_epoc = round(Fs*wnd_train(1)):round(Fs*wnd_train(2)); 
    train_epoc = train_epoc(1:end-mod(length(train_epoc),wnd_test)); % make it divisible
    for k = 1:2
        idxs = find(mrk==k);
        nclass{k} = numel(idxs); % # trials from class k
        for trial = 1:nclass{k}
            tmp_data = eTrain(idxs(trial)+train_epoc,:);
            train_data{k,trial} = tmp_data;
            tmp_C = tmp_data'*tmp_data; % sample covariance
            C{k}(:,:,trial) = tmp_C./trace(tmp_C);
        end
        avrC{k} = mean(C{k}(:,:,1:nclass{k}),3); % class average covariance
    end
    
%     wnd_test = 0.5*Fs;
%     nstep = wnd_test/2; % allow overlap of window
%     train_epoc = round(Fs*wnd_train(1)):round(Fs*wnd_train(2)); 
%     train_epoc = train_epoc(1:end-mod(length(train_epoc),nstep)); % make it divisible
%     for k = 1:2
%         idxs = find(mrk==k);
%         i = 0;
%         for trial = 1:numel(idxs)
%             tmp_data = eTrain(idxs(trial)+train_epoc,:);
%             for j = 1:nstep:size(tmp_data,1)-wnd_test %%%
%                 i = i+1;
%                 tmp_data_i = tmp_data(j:j+wnd_test-1,:);
%                 train_data{k,i} = tmp_data_i;
%                 tmp_C = tmp_data_i'*tmp_data_i; % sample covariance
%                 C{k}(:,:,i) = tmp_C./trace(tmp_C);
%             end
%         end
%         avrC{k} = mean(C{k},3); % class average covariance
%         nclass{k} = i; % # trials from class k 
%     end

    %% CSP
    [Uvec,Uval] = eig(avrC{1}+avrC{2});
    [Uval,idxs] = sort(diag(Uval),'descend'); % sort eigenvals in descend order
    Uvec = Uvec(:,idxs);
    P = sqrt(inv(diag(Uval)))*Uvec'; % whitening matrix
    s1 = P*avrC{1}*P';
    [Rvec,Rval] = eig(s1);
    [Rval,idxs] = sort(diag(Rval),'descend'); 
    R = Rvec(:,idxs);
    W = R'*P;
    Spfilt = W([1:nof end-nof+1:end],:)';   
    
%     [V,D] = eig(avrC{1},avrC{1}+avrC{2});
%     Spfilt = V(:,[1:nof end-nof+1:end]);  
    
    %% Log-variance feature extraction
    for k = 1:2
        features{k} = [];
        for trial = 1:nclass{k}
            for section = 1:length(train_epoc)/wnd_test
                interval = (section-1)*wnd_test+[1:wnd_test];
                features{k} = vertcat(features{k},log(var(train_data{k,trial}(interval,:)*Spfilt)));
            end
        end
    end

%     for k = 1:2
%         features{k} = [];
%         for trial = 1:nclass{k}
%             features{k} = vertcat(features{k}, log(var(train_data{k,trial}*Spfilt)));
%         end
%     end
    
    %% Scatter plots 
    if seeplot == 1 && nof == 4
        figure; hold on
        scatter(features{1}(:,1),features{1}(:,8),'r');
        scatter(features{2}(:,1),features{2}(:,8),'b');
        figure; hold on
        scatter(features{1}(:,2),features{1}(:,7),'r');
        scatter(features{2}(:,2),features{2}(:,7),'b');
        figure; hold on
        scatter(features{1}(:,3),features{1}(:,6),'r');
        scatter(features{2}(:,3),features{2}(:,6),'b');
    end
    
    %% LDA
    lda_W = ((mean(features{2})-mean(features{1}))/(cov(features{1})+cov(features{2})))';
    lda_B = (mean(features{1})+mean(features{2}))*lda_W/2;
    
end