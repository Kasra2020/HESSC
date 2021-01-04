function [CSmap, Processing_time] = HESSC_fauto(img, num_lev)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements the automatic veriosn of HESSC, in case the number
% of clusters are unknown.
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
img = double(img);
[nr,nc,nb]=size(img);
X=reshape(img,nr*nc,nb);
mask=sum(X,2)>0;
X=X(mask,:);
X=X';
X=double(X);
ns = size(X,2);
%% HESSC auto 
tstart=tic;
bta=0.1;%beta parameter is used to change the number of clusters.
tta=0;
alpha=0.99;
lambda=1e-2;
tau = 0.2; %tau parameter is used as the thershold for the binary clustering algorithm. Employing closer vlaue to zero, results in obtianing more detalied result.
n_min = 3;
U = {'U_H','std',[]};   
r = 100;
w = ones(r,1); % the weight of each partitioning
rep = 10; % the number of ECC runs
maxIter = 40;
minThres = 1e-5;
utilFlag = 0;
t = tree(1:ns);
errC=zeros(1,2^(num_lev+1));
errP=zeros(1,2^(num_lev+1));
 for iter = 1:num_lev
        inds = t.depthfirstiterator;
        for i = 1:length(inds)
            if t.isleaf(inds(i))==1
                idx = t.get(inds(i));
                if length(idx)>=n_min
                    tmpt=[];
                    for o = 1:r
                    y = lasso_binary_clustering(X(:,idx),tau, lambda);
                    y(y==0) = 2;
                    tmpt = [tmpt;y];
                    disp(o)
                    end  
                    tmpt=tmpt';
                    [~,pi_index1,~,~,~] = RunECC(tmpt,size(unique(y),2),U,w,rep,maxIter,minThres,utilFlag);
                    if tta>=1
                       TC = X(:,idx);
                       tmpC = TC*TC';
                       [eigvector,eigva] = eig(tmpC);
                       [eigvat,~] = sort(eigva,'descend');
                       eigvat = diag(eigvat);
                       E = cumsum(eigvat);
                       E = E./E(end);
                       d = 0;
                       id = [];
                       for k = 1:nb
                           if E(k)<=alpha
                               id = [id,k];
                           else
                               break;        
                           end
                       end
                       if size(id,2)~=0
                           eignV = eigvector(:,1:size(id,2));
                           newX = zeros(size(TP,2),size(id,2));
                       else
                           eignV = eigvector;
                           newX = zeros(size(TP));
                       end
                       newX = TC'*eignV;
                       tmpDC = eignV*eignV';
                       newX = newX';
                       errC(inds(i)) = sum(abs(newX-tmpDC*newX).^2)/sum(abs(newX).^2);
                       if abs(errP(t.getparent(inds(i)))-errC(inds(i)))/abs(errP(t.getparent(inds(i))))>=bta
                        t = t.addnode(inds(i),idx(pi_index1==1));
                        t = t.addnode(inds(i),idx(pi_index1==2));
                        errP(inds(i)) = errC(inds(i));
                       end
                        tta = tta+1;
                    else
                        t = t.addnode(inds(i),idx(pi_index1==1));
                        t = t.addnode(inds(i),idx(pi_index1==2));
                        TP = X(:,idx);
                        tmpP = TP*TP';
                        [VP,~] = eig(tmpP);
                        tmpDP = VP*VP';
                        errP(inds(i)) = sum(abs(TP-tmpDP*TP).^2)/sum(abs(TP).^2);
                        TP = X(:,idx);
                        tmpP = TP*TP';
                        [eigvector,eigva] = eig(tmpP);
%                         [eigvat,I] = sort(eigva,'descend');
                        [eigvat,~] = sort(eigva,'descend');
                        eigvat = diag(eigvat);
                        EP = cumsum(eigvat);
                        EP = EP./EP(end);
                        d=0;
                        id=[];
                        for k = 1:nb
                            if EP(k,1)<=alpha
                                id = [id,k];
                            else
                                break;        
                            end
                        end
                        if size(id,2)~=0
                            eignV = eigvector(:,1:size(id,2));
                            newX = zeros(size(TP,2),size(id,2));
                        else
                            eignV = eigvector;
                            newX = zeros(size(TP));
                        end
                        newX = TP'*eignV;
                        tmpDP = eignV*eignV';
                        newX = newX';
                        errP(inds(i)) = sum(abs(newX-tmpDP*newX).^2)/sum(abs(newX).^2);
                     tta=tta+1;
                    end                
                end
            end
        end
         disp(length(inds))
end

Segs = ones(1,ns);
inds = t.depthfirstiterator;
k = 1;
for i = 1:length(inds)
    if t.isleaf(inds(i))==1
        k = k+1;
        idx = t.get(inds(i));
        Segs(idx) = k;
    end
end
Processing_time = toc(tstart);
CSmap = zeros(1,nr*nc);
CSmap(1,mask) = Segs;
CSmap = reshape(CSmap,nr,nc);

%% Visualization
figure;imagesc(CSmap);axis off;axis image;
