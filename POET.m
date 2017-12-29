% POET.m carries out POET computational routine, introduced in
%   
% Fan, J, Liao, Y., and Mincheva, M. (2013), 'Large covariance estimation by thresholding
% principal orthogonal complements'. Journal of the Royal Statistical Society: Series B
% (Statistical Methodology), 75(4):603–680.
%
% The INPUT arguments are: 
% 'z': an input data matrix nxp (n objects p variables).
% 'r_thr': the input low rank.
% 'rho': a vector of constants for sparsity threshold selection.
% 'CV_ind': defaults to 1. It activates the cross-validation procedure for
% threshold selection.
% 'H': the number of folds for cross-validation.
% 'ad': defaults to 1, in which case it performs adaptive thresholding 
% (Fan et al., 2013). If set to 0, constant thresholding is performed.
% 'hard': defaults to 1, in case which it performs hard thresholding. 
% It should be set to 0 if soft thresholding is desired.
% 'th_ind': 1 if the theoretical parameters are known, 0 otherwise.
% 'A': the true sparse component, if known.
% 'B': the true low rank component, if known.
% 'Sigma': the true covariance matrix, if known.
%
% The OUTPUT arguments are:
% 'L': the estimated low rank component.
% 'S': the estimated sparse component.
% 'Sigma': the estimated covariance matrix.
% 'rho_opt': the optimal cross-validtation constant for sparsity threshold selection.
% 'non-zeros percentage': the estimated percentage of non-zeros.
% 'latent variance percentage': the estimated percentage of latent variance.
% 'residual covariance percentage': the estimated percentage of residual
% covariance.

function[Out_POET]=POET(z,r_thr,const,CV_ind,H,ad,hard,th_ind,A,B,Sigma)

C=cov(z);
M=C;

n=size(z,1);
p=size(z,2);

ad=1;
hard=1;

if th_ind==1
    r=rank(B)
else
    r=r_thr;
end

CV_ind=1;

if CV_ind==1
for h=1:H
N_max=10000;

eig(M);
[U_C,D_C]=svds(M,r);

clear epsFan SSS vvv vv SS

data=z;
[U_zz D_zz]=(eigs(data'*data,r));
zFan=data*U_zz*U_zz';%1/sqrt(n)*z*princomp(z);1/sqrt(n)*
epsFan=data-zFan;
B_C(:,:,h)=U_C*D_C*U_C';

	%%TUNING

n_training=round(n*(1-inv(log(n))));
n_valid=n-n_training;

VV(1:n_training,h)=randperm(n,n_training)';

epsFan=epsFan';

    s_C=epsFan(:,VV(:,h));
    SS(:,:,h)=cov(s_C');

vv=1:n;
vvv=zeros(n_valid,p);
vvv(:,h)=setdiff(vv',VV(:,h));


ss_C=epsFan(:,vvv(:,h));
SSS(:,:,h)=cov(ss_C');

S=SS(:,:,h);
for t1=1:length(const)
for i=1:(p-1)
    for j=(i+1):p
            more_C=s_C;
        if ad==1
        sigmahat(i,j)=mean(more_C(i,:)*more_C(j,:)');
        for t=1:n_training
        thetapre(t)=(more_C(i,t)*more_C(j,t)-sigmahat(i,j))^2;
        end
        theta(i,j)=mean(thetapre);
        ad_thr(i,j)=const(t1)*theta(i,j);
        else
        ad_thr=const(t1);
        end
        if hard==1
            if S(i,j)>ad_thr(i,j)
               S_C_Thr(i,j)=sign(S(i,j))*abs(S(i,j)); 
            else 
               S_C_Thr(i,j)=0;
            end
        else
        S_C_Thr(i,j)=sign(S(i,j))*max(abs(S(i,j))-ad_thr(i,j),0);
        end
    end;
end;

for i=1:p
    S_C_Thr(i,i)=S(i,i);
end; 
for i=2:p
    for j=1:(i-1)
        S_C_Thr(i,j)=S_C_Thr(j,i);
    end;
end;
S_C(:,:,t1,h)=S_C_Thr;


linf_C_s(t1,h)=norm(S_C(:,:,t1,h),Inf);
l2_C_s(t1,h)=norm(B_C(:,:,h));
Sigma_C(:,:,t1,h)=B_C(:,:,h)+S_C(:,:,t1,h);

lsum=sum(diag(B_C(:,:,h)));
ssum=sum(diag(S_C(:,:,t1,h)));
diagtot_C(t1,h)=lsum+ssum;
rappvar_C(t1,h)=lsum/diagtot_C(t1,h);
v=0;
for i=1:(p-1)
    for j=(i+1):p
        if S_C(i,j,t1)~=0
            v=v+1;
        end;
    end;
end;
v;
numvar=p*(p-1)/2;
nz_C(t1,h)=v;

fac=0;
for i=1:p
for j=(i+1):p
fac=fac+abs(B_C(i,j,h));
end;
end;

plus=0;
for i=1:p
for j=(i+1):p
plus=plus+abs(S_C(i,j,t1,h));
end;
end;

tot_C(t1,h)=plus+fac;
rappcorr_C(t1,h)=plus/tot_C(t1,h);

OB_pond_C(t1,h)=max(linf_C_s(t1)/((1-rappvar_C(t1,h))*trace(M)),l2_C_s(t1,h)/(trace(M)*rappvar_C(t1,h)/r));

TL_CV_C(t1,h)=norm(S_C(:,:,t1,h)-SSS(:,:,h),'fro')^2;

defSpS_C(t1,h)=sum(eig((S_C(:,:,t1,h)))<0);
defSpSigma_C(t1,h)=sum(eig((Sigma_C(:,:,t1,h)))<0);
end;
H
%it
end;
end
    
for t1=1:length(const)
crit_C(t1)=mean(TL_CV_C(t1,:));
defSpSigmamean_C(t1)=mean(defSpSigma_C(t1,:));
defSpmean_C(t1)=mean(defSpS_C(t1,:));
nzmean_C(t1)=mean(nz_C(t1,:));
rappcorrmean_C(t1)=mean(rappcorr_C(t1,:));
end;
%end

crit_C=crit_C(1:length(const))
defSpmean_C=defSpmean_C(1:length(const))
if CV_ind==1
[critminC1 ind_1]=min(crit_C(defSpmean_C==0));
[critminminC ind_2]=min(min(crit_C(defSpmean_C==0)));
critminminC
nd=find(crit_C==critminminC)
ind=nd(1)
end;
rhoFan2=const(ind)
ind
Ind=ind;

[U_C,D_C]=svds(M,r);
    B_C=U_C*D_C*U_C';
    S=M-B_C;
    for t1=1:length(rhoFan2)
for i=2:p-1
for j=(i+1):p
        if ad==1
        sigmahat(i,j)=mean(epsFan(i,:)*epsFan(j,:)');
        for t=1:n_training
        thetapre(t)=(epsFan(i,t)*epsFan(j,t)-sigmahat(i,j))^2;
        end
        theta(i,j)=mean(thetapre);
        ad_thr(i,j)=rhoFan2(t1)*theta(i,j);
        else
        ad_thr(i,j)=rhoFan2(t1);
        end
        if hard==1
            if S(i,j)>ad_thr(i,j)
               S_C_Thr(i,j)=sign(S(i,j))*abs(S(i,j)); 
            else 
               S_C_Thr(i,j)=0;
            end
        else
        S_C_Thr(i,j)=sign(S(i,j))*max(abs(S(i,j))-ad_thr(i,j),0);
        end
end;
end;
        
        %%diagonal thresholding
        for i=1:p
            S_C_Thr(i,i)=S(i,i);
        end;
        for i=2:p
            for j=1:(i-1)
                S_C_Thr(i,j)=S_C_Thr(j,i);
            end;
        end;
        S_C(:,:,t1)=S_C_Thr;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if th_ind==1
        Loss_C(t1)=norm(S_C(:,:,t1)-A,'fro')+norm(B_C-B,'fro');
        linf_C(t1)=norm(S_C(:,:,t1)-A,Inf);
        l2_C(t1)=norm(B_C-B);
        end
        linf_C_s(t1)=norm(S_C(:,:,t1),Inf);
        l2_C_s(t1)=norm(B_C);
        OB_C(t1)=max(linf_C(t1),l2_C(t1));
        
        Sigma_C(:,:,t1)=B_C+S_C(:,:,t1);
        Sigma_Inv_C(:,:,t1)=inv(Sigma_C(:,:,t1));
        
        TL_C(t1)=norm(Sigma_C(:,:,t1)-Sigma,'fro');
        TL_C_s(t1)=norm(Sigma_C(:,:,t1)-M,'fro');
        
        cond_C(t1)=cond(Sigma_C(:,:,t1));
        cond_S_C(t1)=cond(S_C(:,:,t1));
        
        lsum=sum(diag(B_C));
        ssum=sum(diag(S_C(:,:,t1)));
        diagtot_C(t1)=lsum+ssum;
        rappvar_C(t1)=lsum/diagtot_C(t1);
        RelErr_C(t1)=power(p,-1/2)*norm(mpower(inv(Sigma),1/2).*Sigma_C(:,:,t1).*mpower(inv(Sigma),1/2)-eye(p),'fro');
        
        v=0;
        for i=1:(p-1)
            for j=(i+1):p
                if S_C(i,j,t1)~=0
                    v=v+1;
                end;
            end;
        end;
        v;
        %s;

    s=0;
    for i=1:(p-1)
        for j=(i+1):p
            if A(i,j)~=0
                s=s+1;
            end;
        end;
    end;

        numvar=p*(p-1)/2;

        nz_C(t1)=v;
        dist_C(t1)=1-v/s;
        
        spC_11=0;spC_01=0;spC_10=0;spC_00=0;
        sppos_11=0;sppos_01=0;sppos_10=0;sppos_00=0;
        sp0posC=0;sp0negC=0;spnullposC=0;spnullnegC=0;
        for i=1:(p-1)
            for j=(i+1):p
                if A(i,j)~=0 && S_C(i,j,t1)~=0
                    spC_11=spC_11+1;
                end;
                if A(i,j)~=0 && S_C(i,j,t1)==0
                    spC_10=spC_10+1;
                end;
                if A(i,j)==0 && S_C(i,j,t1)~=0
                    spC_01=spC_01+1;
                end;
                if A(i,j)==0 && S_C(i,j,t1)==0
                    spC_00=spC_00+1;
                end;
                if A(i,j)>0 && S_C(i,j,t1)>0
                    sppos_11=sppos_11+1;
                end;
                if A(i,j)<0 && S_C(i,j,t1)<0
                    sppos_00=sppos_00+1;
                end;
                if A(i,j)>0 && S_C(i,j,t1)<0
                    sppos_10=sppos_10+1;
                end;
                if A(i,j)<0 && S_C(i,j,t1)>0
                    sppos_01=sppos_01+1;
                end;
                if A(i,j)>0 && S_C(i,j,t1)==0
                    sp0posC=sp0posC+1;
                end;
                if A(i,j)<0 && S_C(i,j,t1)==0
                    sp0negC=sp0negC+1;
                end;
                if A(i,j)==0 && S_C(i,j,t1)>0
                    spnullposC=spnullposC+1;
                end;
                if A(i,j)==0 && S_C(i,j,t1)<0
                    spnullnegC=spnullnegC+1;
                end;
            end;
        end;
        vp_C(t1)=spC_11;
        fp_C(t1)=spC_01;
        fn_C(t1)=spC_10;
        vn_C(t1)=spC_00;
        vpos_C(t1)=vp_C(t1)+fn_C(t1);
        vneg_C(t1)=vn_C(t1)+fp_C(t1);
        vvpos_C(t1)=sppos_11;
        ffpos_C(t1)=sppos_01;
        ffneg_C(t1)=sppos_10;
        vvneg_C(t1)=sppos_00;
        nnpos_C(t1)=sp0posC;
        nnneg_C(t1)=sp0negC;
        posnn_C(t1)=spnullposC;
        negnn_C(t1)=spnullnegC;
        totpos_C(t1)=vvpos_C(t1)+ffneg_C(t1)+nnpos_C(t1);
        totneg_C(t1)=vvneg_C(t1)+ffpos_C(t1)+nnneg_C(t1);
        tot_C(t1)=spC_11+spC_01+spC_10+spC_00;
        errsp_C(t1)=fn_C(t1)+fp_C(t1);
        errspplus_C(t1)=ffneg_C(t1)+ffpos_C(t1);
        err_C(t1)=errsp_C(t1)/numvar;
        errplus_C(t1)=errspplus_C(t1)/(s);
        %vp/(vp+fn)
        %fn_rate(t1,t2)=fn(t1,t2)/vneg(t1,t2);
        %fp_rate(t1,t2)=fp(t1,t2)/vpos(t1,t2);
        sens_C(t1)=vp_C(t1)/(vpos_C(t1));
        spec_C(t1)=vn_C(t1)/(vneg_C(t1));
        senspos_C(t1)=vvpos_C(t1)/(totpos_C(t1));
        posnegrate_C(t1)=ffneg_C(t1)/(totpos_C(t1));
        posnnrate_C(t1)=nnpos_C(t1)/(totpos_C(t1));
        specpos_C(t1)=vvneg_C(t1)/(totneg_C(t1));
        negposrate_C(t1)=ffpos_C(t1)/totneg_C(t1);
        negnnrate_C(t1)=nnneg_C(t1)/totneg_C(t1);
        possens_C(t1)=posnn_C(t1)/vneg_C(t1);
        negsens_C(t1)=negnn_C(t1)/vneg_C(t1);
        v;
        s;
        
                
    posrate_C(t1)=posnegrate_C(t1)+posnnrate_C(t1);
    negrate_C(t1)=negposrate_C(t1)+negnnrate_C(t1);
    nnrate_C(t1)=possens_C(t1)+negsens_C(t1);
    totrate_C(t1)=(totpos_C(t1)*posrate_C(t1)+totneg_C(t1)*negrate_C(t1)+nnrate_C(t1)*vneg_C(t1))/numvar;
        
        
        fac=0;
        for i=1:p
            for j=(i+1):p
                fac=fac+abs(B_C(i,j));
            end;
        end;
        plus=0;
        for i=1:p
            for j=(i+1):p
                plus=plus+abs(S_C(i,j,t1));
            end;
        end;
        tot_C(t1)=plus+fac;
        rappcorr_C(t1)=plus/tot_C(t1);
        
        eigS_C(:,t1)=svds(S_C_Thr,p);
        errA_C(t1)=norm(eigS_C(:,t1)-svds(A,p),'fro');
        eigL_C=svds(B,r);
        errB_C(t1)=norm(svds(M,r)-svds(B,r),'fro');
        eigSigma_C(:,t1)=svds(Sigma_C(:,:,t1),p);
        errSigma_C(t1)=norm(eigSigma_C(t1)-svds(Sigma,p));
        
defSpS_C(t1)=sum(eig((S_C(:,:,t1)))<0);
defSpSigma_C(t1)=sum(eig((Sigma_C(:,:,t1)))<0);

    end;

if th_ind==1
%rank_it=r_hat;
alpha_it_C=rappvar_C(1,1);
rho_it_C=rappcorr_C(1,1);
nz_it_C=nz_C(1,1);
TL_it_C=TL_C(1,1);
TL_s_it_C=TL_C_s(1,1);
Loss_it_C=Loss_C(1,1);
linf_s_it_C=norm(S_C(:,:,1,1),Inf);
l2_s_it_C=norm(B_C);
Loss_it_C=norm(S_C(:,:,1,1)-A,'fro')+norm(B_C-B,'fro');
Loss_A_C=norm(S_C(:,:,1,1)-A,'fro');
Loss_B_C=norm(B_C-B,'fro');
eig_it_C=norm(eig(Sigma_C(:,:,1))-eig(Sigma));
eig_it_A_C=norm(eig(S_C(:,:,1))-eig(A));
eig_it_B_C=norm(eig(B_C)-eig(B));

end

if defSpS_C(1)==0

Out_POET=cell(1,8);

Out_POET={'L',B_C,'S',S_C(:,:,1),'Sigma',Sigma_C(:,:,1),'rho_opt',rhoFan2,'non-zeros percentage',nz_C(ind)/numvar,'latent variance percentage',rappvar_C(ind),'residual covariance percentage',rappcorr_C(ind)};

else

print("The sparse estimate is not positive definite")

end

end

