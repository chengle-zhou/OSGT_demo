clear all;
close all;
clc

%% load HSI，Lidar dataset

addpath My_SVM
% addpath Efficient Minimization Method for a Generalized Total Variation Functional
load 2013_IEEE_GRSS_DF_Contest_LiDAR.mat
img_Lidar = double(imggt);
load HustonU_gt
load HustonU_IM

%% estimate the size of the input image 
tim2 = hustonu_gt;
tim1 = double(tim2);
[xi,yi] = find(tim1 == 0);
xisize = size(xi);

img = hustonu;
% estimate the size of the input image
[no_lines, no_rows, no_bands] = size(img); 

tim3 = reshape(tim2,1,no_lines*no_rows);
tim3 = tim3';
% vectorization(向量化)
img = ToVector(hustonu);
img = img';

% %% 超像素分割
% 
img_reduction = compute_mapping(img','PCA',3);


img_pca = img_reduction(:,1:3);
img_pca=reshape(img_pca,[no_lines, no_rows, 3]);
max_img_pca = max(max(max(img_pca)));
img_pca = img_pca ./(max_img_pca);
img_pca=ToVector(img_pca);
max_img_Lidar = max(max(img_Lidar));
img_Lidar = img_Lidar ./(max_img_Lidar);
% [img_pca] = scale_new(img_pca);
% [img_Lidar] = scale_new(img_Lidar);
% pca1=img_pca(:,1);
% pca1=reshape(pca1,[no_lines, no_rows, 1]);
% imtool(pca1);

superpixel_data = reshape(img_pca,[no_lines, no_rows, 3]);
superpixel_data=mat2gray(superpixel_data);
superpixel_data=im2uint8(superpixel_data);

superpixel_num= 500;
lambda_prime = 0.8;  sigma = 10;  conn8 = 1;
SuperLabels = mex_ers(double(superpixel_data),superpixel_num,lambda_prime,sigma,conn8);

%% Superpixel Measure Lidar and PCA1
img_PCA = img_pca(:,1);
img_pca1=reshape(img_PCA,[no_lines,no_rows,1]);
img_PCA1=img_PCA';
img_LIDAR=reshape(img_Lidar,[no_rows*no_lines,1]);
img_LIDAR=img_LIDAR';

img_fus = zeros(no_rows*no_lines, 1);
for i=0:superpixel_num-1  %对每个超像素块进行操作
    [SupX_index,SupY_index] = find(SuperLabels==i);
    meas_lables=SuperLabels(min(SupX_index):max(SupX_index),min(SupY_index):max(SupY_index));%将超像素块填充为一个规则矩阵
    meas_PCA=img_pca1(min(SupX_index):max(SupX_index),min(SupY_index):max(SupY_index));%规则矩阵对应的pca data
    meas_PCA=meas_PCA(:);
    number_pixels = size(meas_PCA);
    [xi,yi]=size(meas_lables);
    
    Temp_matri = ones(number_pixels);
    Temp_matri=Temp_matri(:);
    Temp_Matri = ones(number_pixels);
    Temp_Matri=Temp_Matri(:);
    %块内均值
    sup_ind=find(SuperLabels==i);       %key point : position of total HSI 超像素块内每个像素对应在HSI中的位置
    m1=img_PCA1(:,sup_ind);
    supPC_means=mean(m1,2);
    
    position_member = [];
    for j=1:number_pixels  %规则矩阵中属于超像素的像素赋pca data,否则赋均值
        if (meas_lables(j)==i)
            Temp_matri(j) = meas_PCA(j);    % position of location
            position_member = [position_member j];  % the position of superpixel in patched matrix  规则矩阵中属于超像素的位置
        else
            Temp_matri(j) =  supPC_means;
        end
    end
    I=reshape(Temp_matri,[xi,yi,1]);%I是PCA_1超像素块
    %注意，Temp每一次处理一个超像素块的时候会被处理掉
    
    meas_Lidar=img_Lidar(min(SupX_index):max(SupX_index),min(SupY_index):max(SupY_index));
    meas_Lidar=meas_Lidar(:);
    n=img_LIDAR(:,sup_ind);
    supLi_means=mean(n,2);
    
    for h=1:number_pixels
        if (meas_lables(h)==i)
            Temp_Matri(h) = meas_Lidar(h);
        else
            Temp_Matri(h) =  supLi_means;
        end
    end
    V=reshape(Temp_Matri,[xi,yi,1]);
    
    %GTF
    nmpdef;
    pars_irn = irntvInputPars('l1tv');
    
    pars_irn.adapt_epsR   = 1;
    pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
    pars_irn.adapt_epsF   = 1;
    pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff
    pars_irn.pcgtol_ini = 1e-4;
    pars_irn.loops      = 5;
    pars_irn.U0         = I-V;
    pars_irn.variant       = NMP_TV_SUBSTITUTION;
    pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
    pars_irn.pcgtol_ini    = 1e-2;
    pars_irn.adaptPCGtol   = 1;
    
    U = irntv(I-V, {},5, pars_irn);
    X=U+V;
    
    %将超像素块复原为HSI尺寸大小 349*1905
    X_vector = reshape(X, [1, xi*yi]);
    [w, sup_num_matrix] = size(m1); % the number of superpixel pixels in matrix
    %     for z=1:sup_num_matrix
    %         img_fus(sup_ind(z)) = X_vector(position_member(z));
    %     end
    img_fus(sup_ind) = X_vector(position_member);
end

%% ERS guide Lidar and PCA1 fusion image
% Generate mean feature map
img_fus=reshape(img_fus,[no_lines,no_rows,1]);

%% guided filter for each band
% guidance_image = mean_matix;
guidance_image = img_fus;

for i = 1:no_bands
 input_image = hustonu(:,:,i);
r =2; 
eps =0.2^2; 

guided_band_feature =guidedfilter(guidance_image, input_image, r, eps);
guided_feature(:,:,i) = guided_band_feature;
end

final_feature = reshape(guided_feature,no_lines*no_rows, no_bands);

img = final_feature';
%% construct training and test datasets
% set the number of training sample
 train_num = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]; % the paper use abuout 1% 
indexes = train_random_select(GroundT(2,:),train_num); % based on 24 for each class
% get the training-test indexes
train_SL = GroundT(:,indexes);
test_SL = GroundT;
test_SL(:,indexes) = [];

% get the training-test samples and labels
train_samples = img(:,train_SL(1,:))';
train_labels = train_SL(2,:)';
test_samples = img(:,test_SL(1,:))';
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

SVMresult = svmpredict(ones(no_lines*no_rows,1),img,model); 

% Evaluation the performance of the SVM
GroudTest = double(test_labels(:,1));
SVMResultTest = SVMresult(test_SL(1,:),:);
[OA,AA,kappa,CA] = confusion(GroudTest,SVMResultTest);

