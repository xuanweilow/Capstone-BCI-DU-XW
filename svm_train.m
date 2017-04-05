%  https://au.mathworks.com/help/stats/support-vector-machines-for-binary-classification.html
function svm_train( data1,data2,class1,class2,fig )
    
%     length_min = min([length(data1) length(data2)]);
%     data_in = [data1(length(data1)-length_min+1:end)

    data_in = [data1;data2];
    class = [class1.*ones(length(data1),1);class2.*ones(length(data2),1)];
    svmModel = fitcsvm(data_in,class,'KernelFunction','linear','Standardize',true,'ClassNames',[-1,1]);
    
    % Predict scores over the grid
    d = 1;
    [x1Grid,x2Grid] = meshgrid(min(data_in(:,1)):d:max(data_in(:,1)),...
                      min(data_in(:,2)):d:max(data_in(:,2)));
    xGrid = [x1Grid(:),x2Grid(:)];
    [~,scores] = predict(svmModel,xGrid);
   
    % Plot the data and the decision boundary
    if fig == 1 
        figure;
        h(1) = scatter(data1(:,1),data1(:,2),'b','.'); 
        hold on
        h(2) = scatter(data2(:,1),data2(:,2),'r','.');
        h(3) = plot(data_in(svmModel.IsSupportVector,1),data_in(svmModel.IsSupportVector,2),'ko');
        contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
        legend(h,{sprintf('%d',class1),sprintf('%d',class2),'Support Vectors'});
        axis equal
        hold off
    end
    
end
