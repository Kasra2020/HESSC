function [CSmap, Processing_time] = HESSC_fm(img, n_cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements the manual veriosn of HESSC, in case the numeber of
% clusters are know prior to the clustering procedure.
% HESSC was illustrated in the following Paper.
%
% Rafiezadeh Shahi, K., Khodadadzadeh, M., Tusa, L., Ghamisi, P., Tolosana-Delgado, R., & Gloaguen, R. (2020). Hierarchical Sparse Subspace Clustering (HESSC): An Automatic Approach for Hyperspectral Image Analysis. 
% Remote Sensing, 12(15), 2421.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input
% img: Hyperspectral image (3D matrix).
% num_lev: The depth that HESSC needs to proceed (2^{num_lev} identifies the maximum number of clusters)
% Example: setting num_lev = 3; bta = 0; leads to obtian 8 clusters.
%
% output
% CSmap: Clustering map
%
% Example
% See Demo

%% Data preperation
img=double(img);
lambda=1e-2;
[nr,nc,nb]=size(img);
X=reshape(img,nr*nc,nb);
mask=sum(X,2)>0;
X=X(mask,:);
X = X';
ns=size(X,2);

%% HESSC manual structure
tau = 0.2; % 0<tau<1 The less tau the more details
num_lev = round(log2(n_cluster)); % 2^num_lev equals to the max No.clusters
n_min = 3;
U = {'U_H','std',[]};   
r = 100;
w = ones(r,1); % the weight of each partitioning
rep = 10; % the number of ECC runs
maxIter = 40;
minThres = 1e-5;
utilFlag = 0;
t = tree(1:ns);
indx=[];
tstart=tic;
for iter = 1:num_lev
    inds = t.depthfirstiterator;
    for i = 1:length(inds)
        if t.isleaf(inds(i))==1
            indx=[indx,i];
            idx = t.get(inds(i));
            if length(idx)>= n_min
                tmpt=[];
                for o=1:r
                    y = lasso_binary_clustering(X(:,idx),tau,lambda);
                    y(y==0)=2;
                    tmpt=[tmpt;y];
                    disp(o)
                end
                tmpt=tmpt';
                [pi_sumbest,pi_index1,pi_converge,pi_utility,ty] = RunECC(tmpt,2,U,w,rep,maxIter,minThres,utilFlag);
                if(length(indx)<n_cluster)
                    t = t.addnode(inds(i),idx(pi_index1==1));
                    t = t.addnode(inds(i),idx(pi_index1==2));
                else
                    break;
                end             
            end
        end
    end
    disp(length(inds))
end

Segs = ones(1,ns);
inds = t.depthfirstiterator;
k = 0;
for i = 1:length(inds)
      if t.isleaf(inds(i))==1
        k = k+1;
        idx = t.get(inds(i));
        Segs(idx) = k;
    end
end
Processing_time = toc(tstart);
CSmap=zeros(nr*nc,1);
CSmap(mask)=Segs;
CSmap = reshape(CSmap,nr,nc);

%% Visualization
figure;imagesc(CSmap);axis off;axis image;