% tfkm_rob carries out trimmed factorial k-means algorithm, introduced in
%   
% Farnè, M. and Vouldis, A. T. (2017). "Business models of the banks in the
% Euro Area". No. 2070. ECB Working Paper.
%
% Farnè, M. and Vouldis, A. T. (2021). "Banks' business models in the euro area: 
% a cluster analysis in high dimensions". Ann Oper Res 305, 23–57 (2021).
%  
% This method is a robustified version of the factorial k-means algorithm developed
% in Vichi, M., & Kiers, H. A. (2001). Factorial k-means analysis for two-way data. 
% Computational Statistics & Data Analysis, 37(1), 49-64.
%
% The INPUT arguments are: 
% z: n (objects) times p (variables) data matrix.
% alpha: allowed proportion of outliers. 
% (alpha=0 means that the original factorial k-means is applied).
% r_all: vector of tested number of factors r.
% g_all: vector of tested numbers of clusters c.
% cov_yes=1 if the input is the robust MCD covariance matrix estimate,
% (see Rousseeuw, P. & Driessen, K. A fast algorithm for the minimum 
% covariance determinant estimator. Technometrics. 41, 212-223, 1999)
% 0 if the input is the relative robust correlation matrix estimate.
% r_G: 0 for a rational initialization of loadings, memberships and
% centroids (default), 1 for a random initialization 
% (see Farnè, M. and Vouldis, A. T. (2017,2021) for the details).
% N: number of trials (defaults to 100).

% The OUTPUT arguments are:
% factor_scores: the n times r matrix of factor scores for each object.
% eigenvectors: the p times r matrix of eigenvectors.
% centroids: the c times r matrix of centroids.
% r_yes: optimal rank r obtained by maximizing Hartigan statistics.
% g_yes: optimal number of groups c obtained by maximizing Hartigan statistics.
% alpha: outlier proportion
% member_no: the vector of group membership (outliers=0).
% group_num: clusters' size (excluding outliers).
% out_tfkm: identified outliers' indices.
% n_out_tfkm: number of identified outliers.
% member: the vector of group membership (without excluding outliers).
% group_comp: clusters' size (without excluding outliers).
% model_select: 1 if r and c have been selected by Hartigan's statistic, 
% 0 otherwise
% H: if model_select=1, the matrix of Hartigan statistics for each r-c pair 
% such that r<=c-1.

function[output_all]=tfkm_rob(z,alpha,r_all,g_all,cov_yes,r_G,N)

if~exist('r_G','var')
    r_G=0;
end

if~exist('N','var')
    N=100;
end

%if ~exist('alpha','var')
%    alpha=0.1;
%end

if alpha==0%%~exist('robust','var') &&
    robust=0;
end

if ne(alpha,0)==1%%~exist('robust','var') && 
    robust=1;
end

n=size(z,1);%%number of objects
p=size(z,2);%%number of variables

if cov_yes==1
    %%covariance matrix
    C=robustcov(z);
else
    %%correlation matrix
    C=corrcov(robustcov(z));
end

r_length=length(r_all);
g_length=length(g_all);

if length(g_all)==1
    model_select=0
else
    model_select=1
end

if model_select==1
    
    if r_all(1) > g_all(g_length)-1
        % r must be <= c-1
        output_all=strcat('The rank is larger than the number of clusters minus one')
        return;
    end

% the algorithm starts
for i1=r_all(1):r_all(r_length) %% rank
    for i2=g_all(1):g_all(g_length) %% number of groups

        if i1 <=i2-1
        
% Initialize rank
[U_r, D_r]=svds(C,i1);

% Initialisation: 
R=zeros(p,p);
V=zeros(p,p);
E=zeros(p,p);
v=zeros(p,1);
Random=zeros(p,p,N);  
r=i1;
p=length(C);

clear Random E V

for number=1:N

    K=rand(p)*eye(p); 

    % Here a Gram-Schmidt algorithm is used to derive the orthogonal base of the space defined by K

    E(:,1)=K(:,1);  
        
    for j=2:p
    for i=1:(j-1)
        R(i,j)=dot(K(:,j),E(:,i))/(norm(E(:,i))^2); 
    end
    for i=1:(j-1)
    V(:,i)=R(i,j)*E(:,i);  
    end
   
        for h=1:p
        v(h,1)=sum(V(h,1:(j-1)));  
        end
       E(:,j)=K(:,j)-v;    
    end

    for i=1:p
    E(:,i)=E(:,i)/norm(E(:,i)); 
    end

    Random(:,:,number)=E*U_r;
end

% Initialize n_groups

n_groups=i2;
group=1:n_groups;
RawData=z;
% in case we want to run the clustering in a subset of z, this can be
% easily done here.
n_Raw=size(RawData,1);

%% Random group initializer

%r_G=0;
clear G Random_G centr maxRaw minRaw rangeRaw dist_G t_hotelling_pre F_pre pace

%% The clustering initialisation loop (Step 0)
for number=1:N
   
for i=1:p
maxRaw(i)=max(RawData(:,i));      
minRaw(i)=min(RawData(:,i));
rangeRaw(i)=maxRaw(i)-minRaw(:,i);
end

vvar=var(z');
mmean=mean(z');

if r_G==1
for j=1:n_Raw
    casual(j)=randi(n_groups); 
    G(j,casual(j))=1;
end
end

if r_G==0
    fac_pre=z*Random(:,:,number);
for j=1:n_Raw
     dist_G(:,j)=(fac_pre(j,:)-mean(fac_pre))';
     t_hotelling_pre(:,j)=n_Raw*dist_G(:,j)'*inv(cov(fac_pre))*dist_G(:,j);
end   
    
  q_G=quantile(t_hotelling_pre,2*n_groups-1);

  % The distance of each object from each cluster centroid is calculated
  for j=1:n_Raw
  for m=1:n_groups
  F_pre(j,m)=norm(t_hotelling_pre(j)-q_G(2*m-1));
  end
  end
  
% Here the assignment of objects to clusters is done, based on the T-score
for j=1:n_Raw
for m=1:n_groups
    if min(F_pre(j,:))==F_pre(j,m)
        G(j,m)=1;
    else
        G(j,m)=0;
    end
end
end
end
    
Random_G(:,:,number)=G;

clear G
end

%%

clear YY GG UU_GG F member

n_it=0;
%robust=1;

for number=1:N %% number of trials
G=Random_G(:,:,number);
U_G=Random(:,:,number);
Y_bar=inv(G'*G)*G'*RawData*U_G;
ob_pre=norm(RawData*U_G-G*Y_bar)^2;   
diff=ob_pre;
n_it=0;

% The main clustering loop starts. It runs while diff>0
while diff>0

    %% Step 1: membership allocation

    % For each object j, for each group m, the loss for the object is calculated
    for j=1:n
        
        for m=1:n_groups
            count=zeros(1,n_groups);
            count(m)=1; 
            F(j,m)=norm(z(j,:)*U_G-count*Y_bar)^2; 
        end
        
        % each object is assigned to clusters according to the minimum distance
        for m=1:n_groups
            if min(F(j,:))==F(j,m)
                G(j,m)=1;
            else
                G(j,m)=0;
            end
        end

    end

clear member
for j=1:n
    for m=1:n_groups
    if G(j,m)==1
       member(j)=m;
    end
    end
end

%% Step 2: T-score and outliers' identification

centroids=Y_bar;
fac_scores=z*U_G;
cen_true=G*centroids;

if robust==1
for m=1:n_groups
for j=1:n
    diff_scores(:,j)=fac_scores(j,:)'-cen_true(member(j),:)';
    t_hotelling(j)=n*diff_scores(:,j)'*inv(cov(fac_scores))*diff_scores(:,j);
    t_test(j)=(n-r)/(r*(n-1))*t_hotelling(j);
end
end

% the 100(1-alpha)-th percentile of the T-score is calculated
thr_clust_max=quantile(t_test,1-alpha);

clear flag
bound_max=thr_clust_max;
% objects exceeding the T-score threshold are flagged as outliers
for j=1:n
    if t_test(j)>bound_max 
       flag(j)=1;
    else flag(j)=0;
    end
end

% outliers are removed from the groups
for j=1:n
    if flag(j)==1
       G(j,:)=zeros(1,n_groups);
    end
end
end

%% Step 3: updating latent directions and centroids
if rank(G)==size(G,2)  %% G must be full column rank (no void groups)
    
% extraction of eigenvectors in Step 2 
% see ten Berge, J. M. (1993). Least squares optimization in multivariate
% analysis. Leiden, The Netherlands: DSWO Press, Leiden University.
[U_G, D_G]=svds(z'*(G*inv(G'*G)*G'-eye(n))*z,i1);
Y_bar=inv(G'*G)*G'*z*U_G;

%% Step 4: checking the objective function value
if n_it>1
ob_pre=ob_post;
end
ob_post= norm(z*U_G-G*Y_bar)^2;
diff=ob_post-ob_pre;
diff_pre=diff;
n_it=n_it+1;
if diff<0
YY(:,:,number)=Y_bar;
GG(:,:,number)=G;
UU_GG(:,:,number)=U_G;
clear obob
obob(number)=ob_post;
diffdiff(number)=diff;
N_it(number)=n_it;
end
else diff=-1;           
    % if G is not full column rank, it exits immediately, without saving nothing, 
    % and goes on with the next initialiser
end
end
if diff==0
   obob(number)=0;
   N_it(number)=0;
end

end

% We choose the best solution into the set of N initial starts
if sum(N_it)~=0
gg=find(obob==min(obob(N_it~=0)));
end

n_it=0;
%robust=1;
%thr_clust=1.5;

% The procedure is repeated only for the best solution (gg).
for number=gg(1):gg(1)
G=Random_G(:,:,number);
U_G=Random(:,:,number);
Y_bar=inv(G'*G)*G'*z*U_G;
ob_pre=norm(z*U_G-G*Y_bar)^2;
diff=ob_pre;
n_it=0;
while diff>0
for j=1:n
for m=1:n_groups
    count=zeros(1,n_groups);
    count(m)=1;
    F(j,m)=norm(z(j,:)*U_G-count*Y_bar)^2; 
end

for m=1:n_groups
    if min(F(j,:))==F(j,m)
        G(j,m)=1;
    else
        G(j,m)=0;
    end
end

end

clear member
for j=1:n
    for m=1:n_groups
    if G(j,m)==1
       member(j)=m;
    end
    end
end

centroids=Y_bar;
fac_scores=z*U_G;
cen_true=G*centroids;

if robust==1
clear diff_scores t_hotelling
for m=1:n_groups
for j=1:n
    diff_scores(:,j)=fac_scores(j,:)'-cen_true(member(j),:)';
    t_hotelling(j)=n*diff_scores(:,j)'*inv(cov(fac_scores))*diff_scores(:,j);
    t_test(j)=(n-r)/(r*(n-1))*t_hotelling(j);
end
end

thr_clust_max=quantile(t_test,1-alpha);

clear flag
bound_max=thr_clust_max;
for j=1:n
    if t_test(j)>bound_max
       flag(j)=1;
    else flag(j)=0;
    end
end

for j=1:n
    if flag(j)==1
       G(j,:)=zeros(1,n_groups);
    end
end
end

if rank(G)==size(G,2)
[U_G, D_G]=svds(z'*(G*inv(G'*G)*G'-eye(n))*z,i1);
Y_bar=inv(G'*G)*G'*z*U_G;
if n_it>1
ob_pre=ob_post;
end
ob_post= norm(z*U_G-G*Y_bar)^2;
diff=ob_post-ob_pre;
diff_pre=diff;
n_it=n_it+1;
if diff<0
YY(:,:,number)=Y_bar;
GG(:,:,number)=G;
UU_GG(:,:,number)=U_G;
obob(number)=ob_post;
diffdiff(number)=diff;
end
else diff=-1;
end
end
end

% Here, we calculate within cluster sum of squares.

    clear diff_scores
    fac_good=z*U_G;
    it_yes=0;
    if n_it>1
    for i=1:n
        if sum(G(i,:))>0
           it_yes=it_yes+1;
           diff_scores(i)=(norm(fac_good(i,:)'-Y_bar(member(i),:)'))^2;
        end                        
    end

    W(i1,i2)=(sum(reshape(diff_scores,[],1)));%sqrt

    else W(i1,i2)=0;

    end

i1
i2

clearvars -except gg r_G N robust i1 i2 model_select z r_length g_length r_all g_all alpha C z p n W H;
        
        end
    end
end

for k=g_all(1):g_all(g_length-1)
    for r=r_all(1):r_all(r_length)
    if r <= k-1 % the rank can not exceed the number of cluster minus one
       H(r,k+1)=((W(r,k)/W(r,k+1))-1)*(p-k-1);
       if isnan(H(r,k+1))==1 || isinf(H(r,k+1))==1
          H(r,k+1)=0;  
       end
    else
        H(r,k+1)=0; 
    end
    end
end
H

        if sum(sum(H))==0
        output_all=strcat(['It is impossible to select a partion by Hartigan statistics. Please enlarge the vector of cluster numbers g_all and latent ranks r_all, in a way that for at least one c there are two values of r are s.t. r<=c-1.'])
        return;
        end

[max2 ii1]=max(H)
[max1 i2]=max(max2)
if max1>10
    i1=ii1(i2)
else
    i1=r_all(1)
    i2=r_all(1)+1
end
  
else
    
    i1=r_all(1)
    i2=g_all
    
    if i1>i2-1
        % the rank can not exceed the number of cluster minus one
        output_all=strcat('The rank is larger then the number of clusters minus one')
        return;
    end
end


% OPTIMAL SOLUTION IMPLEMENTATION
% Set the optimal rank
[U_r, D_r]=svds(C,i1);

% Initialisation: 
R=zeros(p,p);
V=zeros(p,p);
E=zeros(p,p);
v=zeros(p,1);
Random=zeros(p,p,N);  
r=i1;
p=length(C);

clear Random E V

for number=1:N

    K=rand(p)*eye(p); 

    % Here a Gram-Schmidt algorithm is used to derive the orthogonal base of the space defined by K

    E(:,1)=K(:,1);  
        
    for j=2:p
    for i=1:(j-1)
        R(i,j)=dot(K(:,j),E(:,i))/(norm(E(:,i))^2); 
    end
    for i=1:(j-1)
    V(:,i)=R(i,j)*E(:,i);  
    end
   
        for h=1:p
        v(h,1)=sum(V(h,1:(j-1)));  
        end
       E(:,j)=K(:,j)-v;    
    end

    for i=1:p
    E(:,i)=E(:,i)/norm(E(:,i)); 
    end

    Random(:,:,number)=E*U_r;
end

% Set the optimal number of groups

n_groups=i2;
group=1:n_groups;
RawData=z;
% in case we want to run the clustering in a subset of z, this can be
% easily done here.
n_Raw=size(RawData,1);

%% Random group initializer

%r_G=0;
clear G Random_G centr maxRaw minRaw rangeRaw dist_G t_hotelling_pre F_pre pace

%% The clustering initialisation loop (Step 0)
for number=1:N
   
for i=1:p
maxRaw(i)=max(RawData(:,i));      
minRaw(i)=min(RawData(:,i));
rangeRaw(i)=maxRaw(i)-minRaw(:,i);
end

vvar=var(z');
mmean=mean(z');

if r_G==1
for j=1:n_Raw
    casual(j)=randi(n_groups); 
    G(j,casual(j))=1;
end
end

if r_G==0
    fac_pre=z*Random(:,:,number);
for j=1:n_Raw
     dist_G(:,j)=(fac_pre(j,:)-mean(fac_pre))';
     t_hotelling_pre(:,j)=n_Raw*dist_G(:,j)'*inv(cov(fac_pre))*dist_G(:,j);
end   
    
  q_G=quantile(t_hotelling_pre,2*n_groups-1);

  % The distance of each object from each cluster centroid is calculated
  for j=1:n_Raw
  for m=1:n_groups
  F_pre(j,m)=norm(t_hotelling_pre(j)-q_G(2*m-1));
  end
  end
  
% Here the assignment of objects to clusters is done, based on the T-score
for j=1:n_Raw
for m=1:n_groups
    if min(F_pre(j,:))==F_pre(j,m)
        G(j,m)=1;
    else
        G(j,m)=0;
    end
end
end
end
    
Random_G(:,:,number)=G;

clear G
end

%%

clear YY GG UU_GG F member

n_it=0;
%robust=1;

for number=1:N %% number of trials
G=Random_G(:,:,number);
U_G=Random(:,:,number);
Y_bar=inv(G'*G)*G'*RawData*U_G;
ob_pre=norm(RawData*U_G-G*Y_bar)^2;   
diff=ob_pre;
n_it=0;

% The main clustering loop starts. It runs while diff>0
while diff>0

    %% Step 1: membership allocation

    % For each object j, for each group m, the loss for the object is calculated
    for j=1:n
        
        for m=1:n_groups
            count=zeros(1,n_groups);
            count(m)=1; 
            F(j,m)=norm(z(j,:)*U_G-count*Y_bar)^2; 
        end
        
        % each object is assigned to clusters according to the minimum distance
        for m=1:n_groups
            if min(F(j,:))==F(j,m)
                G(j,m)=1;
            else
                G(j,m)=0;
            end
        end

    end

clear member
for j=1:n
    for m=1:n_groups
    if G(j,m)==1
       member(j)=m;
    end
    end
end

%% Step 2: T-score and outliers' identification

centroids=Y_bar;
fac_scores=z*U_G;
cen_true=G*centroids;

if robust==1
for m=1:n_groups
for j=1:n
    diff_scores(:,j)=fac_scores(j,:)'-cen_true(member(j),:)';
    t_hotelling(j)=n*diff_scores(:,j)'*inv(cov(fac_scores))*diff_scores(:,j);
    t_test(j)=(n-r)/(r*(n-1))*t_hotelling(j);
end
end

% the 100(1-alpha)-th percentile of the T-score is calculated
thr_clust_max=quantile(t_test,1-alpha);

clear flag
bound_max=thr_clust_max;
% objects exceeding the T-score threshold are flagged as outliers
for j=1:n
    if t_test(j)>bound_max 
       flag(j)=1;
    else flag(j)=0;
    end
end

% outliers are removed from the groups
for j=1:n
    if flag(j)==1
       G(j,:)=zeros(1,n_groups);
    end
end
end

%% Step 3: updating latent directions and centroids
if rank(G)==size(G,2)  %% G must be full column rank (no void groups)
    
% extraction of eigenvectors in Step 2 
% see ten Berge, J. M. (1993). Least squares optimization in multivariate
% analysis. Leiden, The Netherlands: DSWO Press, Leiden University.
[U_G, D_G]=svds(z'*(G*inv(G'*G)*G'-eye(n))*z,i1);
Y_bar=inv(G'*G)*G'*z*U_G;

%% Step 4: checking the objective function value
if n_it>1
ob_pre=ob_post;
end
ob_post= norm(z*U_G-G*Y_bar)^2;
diff=ob_post-ob_pre;
diff_pre=diff;
n_it=n_it+1;
if diff<0
YY(:,:,number)=Y_bar;
GG(:,:,number)=G;
UU_GG(:,:,number)=U_G;
clear obob
obob(number)=ob_post;
diffdiff(number)=diff;
N_it(number)=n_it;
end
else diff=-1;           
    % if G is not full column rank, it exits immediately, without saving nothing, 
    % and goes on with the next initialiser
end
end
if diff==0
   obob(number)=0;
   N_it(number)=0;
end

end

% We choose the best solution into the set of N initial starts
if sum(N_it)~=0
gg=find(obob==min(obob(N_it~=0)));
end

n_it=0;
%robust=1;
%thr_clust=1.5;

% The procedure is repeated only for the best solution (gg).
for number=gg(1):gg(1)
G=Random_G(:,:,number);
U_G=Random(:,:,number);
Y_bar=inv(G'*G)*G'*z*U_G;
ob_pre=norm(z*U_G-G*Y_bar)^2;
diff=ob_pre;
n_it=0;
while diff>0
for j=1:n
for m=1:n_groups
    count=zeros(1,n_groups);
    count(m)=1;
    F(j,m)=norm(z(j,:)*U_G-count*Y_bar)^2; 
end

for m=1:n_groups
    if min(F(j,:))==F(j,m)
        G(j,m)=1;
    else
        G(j,m)=0;
    end
end

end

clear member
for j=1:n
    for m=1:n_groups
    if G(j,m)==1
       member(j)=m;
    end
    end
end

centroids=Y_bar;
Y_bar_opt=centroids;
fac_scores=z*U_G;
cen_true=G*centroids;

if robust==1
clear diff_scores t_hotelling
for m=1:n_groups
for j=1:n
    diff_scores(:,j)=fac_scores(j,:)'-cen_true(member(j),:)';
    t_hotelling(j)=n*diff_scores(:,j)'*inv(cov(fac_scores))*diff_scores(:,j);
    t_test(j)=(n-r)/(r*(n-1))*t_hotelling(j);
end
end

thr_clust_max=quantile(t_test,1-alpha);

clear flag
bound_max=thr_clust_max;
for j=1:n
    if t_test(j)>bound_max
       flag(j)=1;
    else flag(j)=0;
    end
end

for j=1:n
    if flag(j)==1
       G(j,:)=zeros(1,n_groups);
    end
end
end

if rank(G)==size(G,2)
U_G_opt=U_G;
[U_G, D_G]=svds(z'*(G*inv(G'*G)*G'-eye(n))*z,i1);
Y_bar=inv(G'*G)*G'*z*U_G;
if n_it>1
ob_pre=ob_post;
end
ob_post= norm(z*U_G-G*Y_bar)^2;
diff=ob_post-ob_pre;
diff_pre=diff;
n_it=n_it+1;
if diff<0
YY(:,:,number)=Y_bar;
GG(:,:,number)=G;
UU_GG(:,:,number)=U_G;
obob(number)=ob_post;
diffdiff(number)=diff;
end
else diff=-1;
    output_all=strcat('There is at least one void group in the recovered solution.')
    return;
end
end
end

%% OUTPUT

if robust==1
member_no=member;
for j=1:n
    if flag(j)==1
       member_no(j)=0;
    end
end

out_tfkm=find(member_no==0);
n_out_tfkm=length(out_tfkm);

clear group_num % gives number of members of each group, excluding outliers!
for m=1:n_groups
group_num(m)=length(find(member_no==m));
end
group_num;
end

% here the composition also including outliers
clear group_comp
for m=1:n_groups
group_comp(m)=length(find(member==m));
end

if robust==0
    out_tfkm=0;
    n_out_tfkm=0;
    group_num=group_comp;
    member_no=member;
end

   r_yes=i1;
   g_yes=i2;

if model_select==0
output_all={fac_scores,'factor_scores',U_G_opt,'eigenvectors',Y_bar_opt,'centroids',r_yes,'optimal_rank',g_yes,'optimal_cluster_number',alpha,'alpha',member_no,'member_no',group_num,'group_num',out_tfkm,'out_tfkm',n_out_tfkm,'n_out_tfkm',member,'member',group_comp,'group_comp',model_select,'model_select'}
end

if model_select==1
output_all={fac_scores,'factor_scores',U_G_opt,'eigenvectors',Y_bar_opt,'centroids',r_yes,'optimal_rank',g_yes,'optimal_cluster_number',alpha,'alpha',member_no,'member_no',group_num,'group_num',out_tfkm,'out_tfkm',n_out_tfkm,'n_out_tfkm',member,'member',group_comp,'group_comp',model_select,'model_select',H,'Hartigan'}
end