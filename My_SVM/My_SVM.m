function Muti_fea_results = My_SVM(RF_feature,mean_matix,weighted_matix,train_SL,test_SL,rows, cols)
% Author: Chengle Zhou
Muti_fea_results = [];
R_fea = normcol(RF_feature'); % bands * (rows*cols)
M_fea = normcol(mean_matix'); % bands * (rows*cols)
W_fea = normcol(weighted_matix'); % bands * (rows*cols)

fea_sel = {'RF_feature','mean_matix','weighted_matix'};
for i = 1 : length(fea_sel)
    temp_img = eval(fea_sel{i});
    img = temp_img';
    % get the training-test samples and labels
    train_samples = img(:,train_SL(1,:))';
    train_labels = train_SL(2,:)';
    test_labels = test_SL(2,:)';
    
    % Normalize the training set and original image
    [train_samples,M,m] = scale_func(train_samples);
    [img ] = scale_func(img',M,m);
    % Select the paramter for SVM with five-fold cross validation
    [Ccv Gcv cv cv_t] = cross_validation_svm(train_labels,train_samples);
    % Training using a Gaussian RBF kernel
    % give the parameters of the SVM (Thanks Pedram for providing the
    % parameters of the SVM)
    parameter = sprintf('-c %f -g %f -m 500 -t 2 -q',Ccv,Gcv);    
    % Train the SVM
    model = svmtrain(train_labels,train_samples,parameter);
    % SVM Classification
    SVMresult = svmpredict(ones(rows*cols,1),img,model);
    Muti_fea_results = [Muti_fea_results,SVMresult];
%     % Evaluation the performance of the SVM
%     GroudTest = double(test_labels(:,1));
%     SVMResultTest = SVMresult(test_SL(1,:),:);
%     % [SVMOA,SVMAA,SVMkappa,SVMCA]=confusion(GroudTest,SVMResultTest)
%     [OA,AA,kappa,CA] = confusion(GroudTest,SVMResultTest);

end
