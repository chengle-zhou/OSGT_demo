 function [data M m] =scale_func(data,M,m)
%
% function data =rescale(data)
                          %rescale（A）将数组的条目缩放到区间 [0,1]
%
% This function rescale the input data between -1 and 1
%
% INPUT
%
% data: the data to rescale
% max: the maximum value of the ouput data
% min: the minimum value of the output data
% 
% OUTPUT
%
% data: the rescaled data
[Nb_s Nb_b]=size(data);
if nargin==1  %nargin函数输入参数的个数
    M=max(data,[],1);%返中回data每一列的最大值
    m=min(data,[],1);%返中回data每一列的最小值
end

data = 2*(data-repmat(m,Nb_s,1))./(repmat(M-m,Nb_s,1))-1;
%repmat（A,m,n）以A的内容堆叠在（MxN）的矩阵B中，B矩阵的大小由MxN及A矩阵的内容决定