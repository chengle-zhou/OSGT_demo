function[irnPars] = irntvInputPars(dname)
%  
% function [irnPars] = irntv(methodName)
%
% Genrates de defaults parameters to be used in the IRN method, which
% computes the minimum of a generalised TV functional 
%
%           T = || K*U - S ||^p + lambda*|| sqrt( (Dx(U))^2 + (Dy(U))^2 ) ||^q
%               s.t. 0 <= U <= vmax
%  
%           for grayscale / color (vector) images using the IRN [1,2,3,4] algorithm,
%           and p, q in (0 2]
%
%
%           [1] "Efficient Minimization Method for a Generalized Total 
%                Variation Functional"
%               IEEE Transactions on Image Processing, 2009, 18:2(322-332)
%
%           [2] "A Generalized Vector-Valued Total Variation Algorithm"
%               Proceedings of the 2009 International Conference on
%               Image Processing (ICIP 2009), pp. 1309-1312, 2009
%  
%           [3] "A Non-negative Quadratic Programming approach to Minimize  
%               the Generalized Vector-Valued Total Variation Functional"
%               Proceedings of the 2010 European Signal Processing Conference
%                (EUSIPCO 2010), pp. 314-318, 2010
%  
%           [4] "Multiplicative Updates Algorithm To Minimize The Generalized 
%               Total Variation Functional With A Non-Negativity Constraint"
%               Proceedings of the 2010 International Conference on
%               Image Processing (ICIP 2010), pp. 2509-2512, 2010
%  
%           irntv.m is part of NUMIPAD (http://numipad.sf.net). 
%
%           The NUMIPAD library is being developed under U.S. Government
%           contract W-7405-ENG-36 for Los Alamos National Laboratory.
%
%
% Usage:
%       irnPars = irntvInputPars(methodName)
%
% Input:
%       methodName  'l1tv',         l1-TV denoising/deconvolution
%                   'l2tv',         l2-TV denoising/deconvolution
%                   'l1tv_nqp',     l1-TV denoising/deconvolution with 
%                                   a non-negativity variant
%                   'l2tv_nqp'      l2-TV denoising/deconvolution with
%                                   a non-negativity variant
%
% Output:
%       irnPars     Structure with parameters for the IRN method
%       
%       irnPars.p             Data fidelity term norm   数据保真度
%       irnPars.q             Data regularization term norm  正则化项
%       irnPars.loops         Number of iterations  迭代次数
%       irnPars.U0            Initial solution  初始解决方案
%       irnPars.lambda_ini    Regularisation term weighting factor for the 
%                             initial solution (valid for 'l1tv' and/or 'l2tv')初始解的正则项加权因子
%       irnPars.pcgtol_ini    PCG tolerance for the initail solution (valid for初始解决方案的PCG公差
%                             'l1tv' and/or 'l2tv')
%       irnPars.pcgitn        Maximum number of iteration for PCG (valid for      PCG的最大迭代次数
%                             'l1tv' and/or 'l2tv')
%       irnPars.rrs           Relative residual scaling after each iteration每次迭代后的相对残差缩放
%       irnPars.sbstflg       Substitution flag. For p ~= 2 this options 
%                             improve performance. (see [1]). 替代标志。 对于p?= 2，此选项可提高性能。 （请参阅[1]）
%       irnPars.epsf          Data fidelity epsilon (not used if   数据保真度epsilon（如果设置了irnPars.adapt_epsF，则不使用）
%                             irnPars.adapt_epsF is set)
%       irnPars.epsr          Regularization epsilon (not used if   正则化epsilon（如果设置了irnPars.adapt_epsR，则不使用）
%                             irnPars.adapt_epsR is set)
%       irnPars.adapt_epsF    Auto-adapt the data fidelity epsilon. See [1], 
%                             section G.      自动适应数据保真度epsilon。 参见[1]，G节。
%       irnPars.adapt_epsR    Auto-adapt the data regularization epsilon. See [1], 
%                             section G.     自动适应数据正则化epsilon。 参见[1]，G节。
%       irnPars.epsR_cutoff   Cut off value when using the Auto-adapt. 
%                             See [1], section G.    使用自动调整时，请切断值。 参见[1]，G节。
%       irnPars.epsF_cutoff   Cut off value when using the Auto-adapt. 
%                             See [1], section G.
%  
%       irnPars.loops_NQP     Maximun number of iterations for the NQP solver (valid
%                             for 'l1tv_nqp' and/or 'l2tv_nqp')    NQP求解器的最大迭代次数
%       irnPars.vmax_NQP      Value of the variant (default Inf).    变量的值（默认Inf）        
%       irnPars.gamma_NQP     Tolerance parameter to break the inner loop of the NQP 
%                             solver. See section 2.4 of [3] (valid for 'l1tv_nqp' 
%                             and/or 'l2tv_nqp')   公差参数可打破NQP求解器的内部循环。 请参阅[3]的2.4节
%       irnPars.alpha_NQP     Tolerance parameter to break the inner loop of the NQP 
%                             solver. See section 2.4 of [3] (valid for 'l1tv_nqp' 
%                             and/or 'l2tv_nqp')
%                             
%  
% Example:
%       irnPars = irntvInputPars('l1tv')
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov

if(nargin == 0)
  dname='none';
end

nmpdef;

irnPars = struct('p', [], 'q', [], 'epsR', [], 'epsF', [], ...
                 'loops', [], 'U0', [], 'lambda', [], 'lambda_ini', [], ...
                 'pcgtol_ini', [], 'pcgitn', [], 'constraint', [], ...
                 'rrs', [], 'sbstflg', [], ...
                 'weight_scheme', [], ...
                 'adapt_epsR', [], 'adapt_epsF', [], ...
                 'epsR_cutoff', [], 'epsF_cutoff', [], ...
                 'loops_NQP', [], ...
                 'gamma_NQP',[], 'alpha_NQP',[], 'vmax_NQP',[], 'variant', [], ...
                  'M', [], 'speckle_variant', [], 'problem', [], 'adaptPCGtol', []);


switch lower(dname)

  case{'lptv'}
    irnPars.adapt_epsR    = 1;
    irnPars.adapt_epsF    = 1;
    irnPars.problem       = NMP_LPTV;
    irnPars.variant       = NMP_TV_SUBSTITUTION;        

  case{'l1tv'}
    irnPars.p             = 1;
    irnPars.q             = 1;
    irnPars.adapt_epsR    = 1;
    irnPars.adapt_epsF    = 1;
    irnPars.problem       = NMP_L1TV;
    irnPars.variant       = NMP_TV_SUBSTITUTION;        
    irnPars.adaptPCGtol   = 1;        

  case{'l1tv_adapt'}
    irnPars.p             = 1;
    irnPars.q             = 1;
    irnPars.sbstflg       = 1;
    irnPars.adapt_epsR    = 1;
    irnPars.adapt_epsF    = 1;
    irnPars.problem       = NMP_L1TV;
    irnPars.variant       = NMP_TV_LamdaAdapt;

  case{'l1tv_l2tv'}
    irnPars.p             = 1;
    irnPars.q             = 1;
    irnPars.sbstflg       = 0;
    irnPars.adapt_epsR    = 1;
    irnPars.adapt_epsF    = 1;
    irnPars.problem       = NMP_L1L2TV;
    irnPars.variant       = NMP_TV_LamdaAdapt;


  case{'l2tv'}
    irnPars.p             = 2;
    irnPars.q             = 1;
    irnPars.sbstflg       = 0;
    irnPars.adapt_epsR    = 1;
    irnPars.problem       = NMP_L2TV;
    irnPars.variant       = NMP_TV_STANDARD;
    irnPars.adaptPCGtol   = 1;        

  case{'l1tv_nqp'}
    irnPars.p           = 1;
    irnPars.q           = 1;
    irnPars.sbstflg     = 0;
    irnPars.loops_NQP   = 10;
    irnPars.gamma_NQP   = -1;
    irnPars.alpha_NQP   = 1;
    irnPars.vmax_NQP    = Inf;
    irnPars.problem     = NMP_L1TV;
    irnPars.variant     = NMP_TV_NQP;

  case{'l2tv_nqp'}
    irnPars.p           = 2;
    irnPars.q           = 1;
    irnPars.sbstflg     = 0;
    irnPars.loops_NQP   = 10;
    irnPars.gamma_NQP   = -1;
    irnPars.alpha_NQP   = 1;
    irnPars.vmax_NQP    = Inf;
    irnPars.problem     = NMP_L2TV;
    irnPars.variant     = NMP_TV_NQP;


  case{'tv_poisson'}
    irnPars.p           = 2;  % Not used
    irnPars.q           = 1;
    irnPars.sbstflg     = 0;
    irnPars.loops_NQP   = 10;
    irnPars.gamma_NQP   = -1;
    irnPars.alpha_NQP   = 1;
    irnPars.vmax_NQP    = Inf;
    irnPars.problem     = NMP_TV_POISSON;

  case{'tv_poisson_adapt'}
    irnPars.p           = 2;  % Not used
    irnPars.q           = 1;
    irnPars.sbstflg     = 0;
    irnPars.loops_NQP   = 10;
    irnPars.gamma_NQP   = -1;
    irnPars.alpha_NQP   = 1;
    irnPars.vmax_NQP    = Inf;
    irnPars.problem     = NMP_TV_POISSON;
    irnPars.variant     = NMP_TV_LamdaAdapt;

  case{'tv_speckle'}
    irnPars.p           = 2;  % Not used
    irnPars.q           = 1;
    irnPars.sbstflg     = 0;
    irnPars.loops_NQP   = 10;
    irnPars.gamma_NQP   = -1;
    irnPars.alpha_NQP   = 1;
    irnPars.vmax_NQP    = Inf;
    irnPars.problem     = NMP_TV_SPECKLE;
    irnPars.M           = 1;
    irnPars.constraint  = NMP_TV_SPECKLE;
    irnPars.variant     = -1; %FIXME

end


irnPars.weight_scheme = NMP_WEIGHTS_THRESHOLD;

irnPars.loops      = 5;
irnPars.pcgtol_ini = 1e-3;
irnPars.pcgitn     = 500;
irnPars.rrs        = 5;
irnPars.epsF       = 1e-2;
irnPars.epsR       = 1e-4;

irnPars.epsF_cutoff    = 0.05;
irnPars.epsR_cutoff    = 0.01;
