% UNALCE.m carries out UNALCE computation routine, introduced in
%   
% Farnè, M. and Montanari, A. (2017), 'A finite sample estimator 
% for large covariance matrices'.
%
% The INPUT arguments are: 
% C: an input covariance matrix estimator (usually the sample one).
% lambda: a vector of spikiness threshold parameters.
% rho: a vector of sparsity threshold parameters.
% UNALCE: defaults to 1. If set to another value, it performs LOREC
% routine.
% th_ind: 1 if the theoretical parameters are known, 0 otherwise.
% A: the true sparse component, if known.
% B: the true low rank component, if known.
% Sigma: the true covariance matrix, if known.
%
% The OUTPUT arguments are:
% L: the estimated low rank component.
% S: the estimated sparse component.
% Sigma: the estimated covariance matrix.
% lambda_opt: the optimal spikiness threshold selected by MC.
% rho_opt: the optimal sparsity threshold selected by MC.
% rank: the estimated low rank. 
% non-zeros percentage: the estimated percentage of non-zeros.
% latent variance percentage: the estimated percentage of latent variance.
% residual covariance percentage: the estimated percentage of residual
% covariance.

function[Out_UNALCE]=UNALCE(C,lambda,rho,UNALCE,th_ind,A,B,Sigma)

    N_max=1000;
    
    p=size(C,1);

    %r=rank(B);
    
    M=C;
    M_orig=C;
    M_star=C;

    k=1;
    arr=zeros(1,N_max);
    al=zeros(1,N_max);
    criterion=zeros(1,N_max);
    arr(1)=1;
    criterion(1)=1;
    al(1)=1;

    L_Thr=diag(diag(M))/2;
    S_Thr=diag(diag(M))/2;
    E=zeros(p,p);
    Y=L_Thr;
    Z=S_Thr;

%% start loop in lambda and rho

    for t1=1:length(rho) 
    for t2=1:length(lambda)

%% start minimization loop
while k<N_max && abs(criterion(k))>1.0e-02
     
    ARG1=Y-1/2*(Y+Z-M);
    [U,D]=svds(ARG1,rank(ARG1));

	%% SVT
    for i=1:rank(D)
        D_Thr(i,i)=max(D(i,i)-lambda(t2),0);
    end;


    r_Thr=rank(D_Thr);
    rank_Thr(t1,t2)=r_Thr;
    L_pre=L_Thr;
    L_Thr=U(:,1:r_Thr)*D_Thr(1:r_Thr,1:r_Thr)*U(:,1:r_Thr)';

	%% error low rank estimate
    add1=norm(L_Thr-L_pre,'fro')/(1+norm(L_pre,'fro'));


	%% soft-thresholding
    S=Z-1/2*(Y+Z-M);
    S_pre=S_Thr;
    M_pre=L_pre+S_pre;

    for i=1:(p)
    for j=i+1:p
        S_Thr(i,j)=sign(S(i,j))*max(abs(S(i,j))-rho(t1),0);
    end;
    end;


    for i=1:p
       S_Thr(i,i)=S(i,i);
    end; 

    for i=2:p
        for j=1:(i-1)
            S_Thr(i,j)=S_Thr(j,i);
        end;
    end;

    %S_Thr;
    %rank(S_Thr);
	%% error sparse estimate
    add2=norm(S_Thr-S_pre,'fro')/(1+norm(S_pre,'fro'));


    M_star=L_Thr+S_Thr;
    E=M_star-M;

	%% convergence criterion
    k=k+1;
    arr(k)=norm(E)/norm(M);
    criterion(k)=add1+add2;

	%% update L and S	
    al(k)=(1+sqrt(1+4*al(k-1)^2))/2;
    Y=L_Thr+((al(k)-1)/al(k-1))*(L_Thr-L_pre);
    Z=S_Thr+((al(k)-1)/al(k-1))*(S_Thr-S_pre);
    end;
%% end minimization loop

    if UNALCE==1
	%% unshrinkage procedure
    for i=1:r_Thr
    D_Thr(i,i)=D_Thr(i,i)+lambda(t2);
    end;
    L_Thr=U(:,1:r_Thr)*D_Thr(1:r_Thr,1:r_Thr)*U(:,1:r_Thr)';

    for i=1:p
    S_Thr(i,i)=M_star(i,i)-L_Thr(i,i);
    end;
    
    end

	
	%% non-zero detection
    v=0;
    for i=1:(p-1)
        for j=(i+1):p
            if S_Thr(i,j)~=0
                v=v+1;
            end;
        end;
    end;
    
    numvar=p*(p-1)/2;

	%% storage and statistics
    nz(t1,t2)=v;
    Low(:,:,t1,t2)=L_Thr;
    Sparse(:,:,t1,t2)=S_Thr;
    rank_Thr(t1,t2)=r_Thr;
    Sigma_hat(:,:,t1,t2)=L_Thr+S_Thr;
    E_New=M-Sigma_hat(:,:,t1,t2);
    W(:,:,t1,t2)=E_New;
    WW(:,:,t1,t2)=E_New'*E_New;

   
    TL_s(t1,t2)=norm(Sigma_hat(:,:,t1,t2)-C,'fro');
    lambda_d=diag(D_Thr);
    if rank(D_Thr)~=0
    condB(t1,t2)=lambda_d(1)/lambda_d(rank(D_Thr));
    end;
    condA(t1,t2)=cond(S_Thr);
    condSigma_hat(t1,t2)=cond(Sigma_hat(:,:,t1,t2));

    t_L(t1,t2)=trace(Low(:,:,t1,t2));
    t_S(t1,t2)=trace(Sparse(:,:,t1,t2));
    t_W(t1,t2)=trace(W(:,:,t1,t2));
    t_TOT(t1,t2)=trace(Low(:,:,t1,t2))+trace(Sparse(:,:,t1,t2));
    rappvar(t1,t2)=t_L(t1,t2)/t_TOT(t1,t2);

    lfro(t1,t2)=norm(W(:,:,t1,t2),'fro');
    l1(t1,t2)=norm(S_Thr,1);
    lnuc(t1,t2)=sum(diag(D_Thr));

    linf_s(t1,t2)=norm(S_Thr,Inf);
    l2_s(t1,t2)=norm(L_Thr);

    Res(:,:,t1,t2)=W(:,:,t1,t2);
    Arr(t1,t2)=k;

    Add_lfro(t1,t2)=1/2*lfro(t1,t2)^2;
    Add_l1(t1,t2)=rho(t1)*l1(t1,t2);
    Add_lnuc(t1,t2)=lambda(t2)*lnuc(t1,t2);
    ob(t1,t2)=Add_lfro(t1,t2)+rho(t1)*Add_l1(t1,t2)+lambda(t2)*Add_lnuc(t1,t2);

    linf_s(t1,t2)=norm(S_Thr,Inf);
    l2_s(t1,t2)=norm(L_Thr);


    if th_ind==1
    
    s=0;
    for i=1:(p-1)
        for j=(i+1):p
            if A(i,j)~=0
                s=s+1;
            end;
        end;
    end;
        
        
    linf(t1,t2)=norm(S_Thr-A,Inf);
    l2(t1,t2)=norm(L_Thr-B);    
    TL(t1,t2)=norm(Sigma_hat(:,:,t1,t2)-Sigma,'fro');
    dist(t1,t2)=1-v/s;
    Loss(t1,t2)=norm(S_Thr-A,'fro')+norm(L_Thr-B,'fro');
    %rankerr(t1,t2)=r_Thr-r;
	%% sparsity pattern detection
    sp_11=0;sp_01=0;sp_10=0;sp_00=0;
    sppos_11=0;sppos_01=0;sppos_10=0;sppos_00=0;
    sp0pos=0;sp0neg=0;spnullpos=0;spnullneg=0;
    for i=1:(p-1)
        for j=(i+1):p
            if A(i,j)~=0 && S_Thr(i,j)~=0
                sp_11=sp_11+1;
            end;
            if A(i,j)~=0 && S_Thr(i,j)==0
                sp_10=sp_10+1;
            end;
            if A(i,j)==0 && S_Thr(i,j)~=0
                sp_01=sp_01+1;
            end;
             if A(i,j)==0 && S_Thr(i,j)==0
                sp_00=sp_00+1;
            end;
            if A(i,j)>0 && S_Thr(i,j)>0
                sppos_11=sppos_11+1;
            end;
            if A(i,j)<0 && S_Thr(i,j)<0
                sppos_00=sppos_00+1;
            end;
            if A(i,j)>0 && S_Thr(i,j)<0
                sppos_10=sppos_10+1;
            end;
            if A(i,j)<0 && S_Thr(i,j)>0
                sppos_01=sppos_01+1;
            end;
            if A(i,j)>0 && S_Thr(i,j)==0
                sp0pos=sp0pos+1;
            end;
            if A(i,j)<0 && S_Thr(i,j)==0
                sp0neg=sp0neg+1;
            end;
            if A(i,j)==0 && S_Thr(i,j)>0
                spnullpos=spnullpos+1;
            end;
            if A(i,j)==0 && S_Thr(i,j)<0
                spnullneg=spnullneg+1;
            end;
        end;
    end;
    vp(t1,t2)=sp_11;
    fp(t1,t2)=sp_01;
    fn(t1,t2)=sp_10;
    vn(t1,t2)=sp_00;
    vpos(t1,t2)=vp(t1,t2)+fn(t1,t2);
    vneg(t1,t2)=vn(t1,t2)+fp(t1,t2);
    vvpos(t1,t2)=sppos_11;
    ffpos(t1,t2)=sppos_01;
    ffneg(t1,t2)=sppos_10;
    vvneg(t1,t2)=sppos_00;
    nnpos(t1,t2)=sp0pos;
    nnneg(t1,t2)=sp0neg;
    posnn(t1,t2)=spnullpos;
    negnn(t1,t2)=spnullneg;
    totpos(t1,t2)=vvpos(t1,t2)+ffneg(t1,t2)+nnpos(t1,t2);
    totneg(t1,t2)=vvneg(t1,t2)+ffpos(t1,t2)+nnneg(t1,t2);
    tot(t1,t2)=sp_11+sp_01+sp_10+sp_00;
    errsp(t1,t2)=fn(t1,t2)+fp(t1,t2);
    errspplus(t1,t2)=ffneg(t1,t2)+ffpos(t1,t2);
    err(t1,t2)=errsp(t1,t2)/numvar;
    errplus(t1,t2)=errspplus(t1,t2)/(s);
    %vp/(vp+fn)
    %fn_rate(t1,t2)=fn(t1,t2)/vneg(t1,t2);
    %fp_rate(t1,t2)=fp(t1,t2)/vpos(t1,t2);
    sens(t1,t2)=vp(t1,t2)/(vpos(t1,t2));
    spec(t1,t2)=vn(t1,t2)/(vneg(t1,t2));
    senspos(t1,t2)=vvpos(t1,t2)/(totpos(t1,t2));
    posnegrate(t1,t2)=ffneg(t1,t2)/(totpos(t1,t2));
    posnnrate(t1,t2)=nnpos(t1,t2)/(totpos(t1,t2));
    specpos(t1,t2)=vvneg(t1,t2)/(totneg(t1,t2));
    negposrate(t1,t2)=ffpos(t1,t2)/totneg(t1,t2);
    negnnrate(t1,t2)=nnneg(t1,t2)/totneg(t1,t2);
    possens(t1,t2)=posnn(t1,t2)/vneg(t1,t2);
    negsens(t1,t2)=negnn(t1,t2)/vneg(t1,t2);
    %v;
    %s
    end;
    
    maximum_eig(t1,t2)=norm(Sigma_hat(:,:,t1,t2));

    M=C;
    M_orig=C;
    M_star=C;
    k=1;
    
	%% algorithm statistics
    Err(1:length(arr),t1,t2)=arr;
    %arr=zeros(1,N_max);
    arr(1)=1;
    
    for i_Crit=1:length(criterion)
    Crit(i_Crit,t1,t2)=criterion(i_Crit);
    end;
    %criterion=zeros(1,N_max);
    criterion(1)=1;

	%% residual covariance proportion
    fac=0;
    for i=1:p
    for j=(i+1):p
    fac=fac+abs(L_Thr(i,j));
    end;
    end;

    plus=0;
    for i=1:p
    for j=(i+1):p
    plus=plus+abs(S_Thr(i,j));
    end;
    end;

    lsum=sum(diag(L_Thr));
    ssum=sum(diag(S_Thr));
    diagtot(t1,t2)=lsum+ssum;
    %rappvar(t1,t2)=lsum/diagtot(t1,t2);

    tott(t1,t2)=plus+fac;
    rappcorr(t1,t2)=plus/tott(t1,t2);

	%% model selection criterion
    
    scale(t1,t2)=rho(t1)/lambda(t2);   
    OB_pond_yes(t1,t2)=max(linf_s(t1,t2)/(scale(t1,t2)*(1-rappvar(t1,t2))*trace(M)),(rank_Thr(t1,t2)*l2_s(t1,t2))/(trace(M)*rappvar(t1,t2)));
    difflinfyes(t1,t2)=OB_pond_yes(t1,t2)-linf_s(t1,t2)/(scale(t1,t2)*(1-rappvar(t1,t2))*trace(M));
    diffl2yes(t1,t2)=OB_pond_yes(t1,t2)-(rank_Thr(t1,t2)*l2_s(t1,t2))/(trace(M)*rappvar(t1,t2));
   
    end;
    end;
    
    [minmin1 ind_1]=min(OB_pond_yes)
    [minmin2 ind_2]=min(min(OB_pond_yes))
    fin1=ind_1(ind_2)
    fin2=ind_2    

Out_UNALCE=cell(1,18);

Out_UNALCE={'L',Low(:,:,fin1,fin2),'S',Sparse(:,:,fin1,fin2),'Sigma',Sigma_hat(:,:,fin1,fin2),'lambda_opt',lambda(fin2),'rho_opt',rho(fin1),'rank',r_Thr,'non-zeros percentage',nz(fin1,fin2)/numvar,'latent variance percentage',rappvar(fin1,fin2),'residual covariance percentage',rappcorr(fin1,fin2)};
end