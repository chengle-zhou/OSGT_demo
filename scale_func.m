 function [data M m] =scale_func(data,M,m)
%
% function data =rescale(data)
                          %rescale��A�����������Ŀ���ŵ����� [0,1]
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
if nargin==1  %nargin������������ĸ���
    M=max(data,[],1);%���л�dataÿһ�е����ֵ
    m=min(data,[],1);%���л�dataÿһ�е���Сֵ
end

data = 2*(data-repmat(m,Nb_s,1))./(repmat(M-m,Nb_s,1))-1;
%repmat��A,m,n����A�����ݶѵ��ڣ�MxN���ľ���B�У�B����Ĵ�С��MxN��A��������ݾ���