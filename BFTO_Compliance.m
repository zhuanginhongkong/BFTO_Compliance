function BFTO_Compliance(lelx, lely, vol_con)
%% PARAMETERS DEFINITION
scale = 1; h = 1;
lelx = 160*scale; lely =100*scale;
BDY = [0, 0; lelx, lely] ;
vol_con = 0.5;
xmin = 1.e-9; dbeta = 1;
E0 = 10; nu = 0.3;
D = E0/(1-nu^2)*[1  nu  0;  nu  1  0;  0  0  (1-nu)/2] ;
rmin =  max(1.2*h, scale * 10 * vol_con); % FILTER RADIUS
Meshdiff = 15*h; d1 = 0.5*h; d2 = 1.0*h;
%% INITIAL MESHING
[xn, yn] = meshgrid(0: h :lelx, 0: h: lely) ;
[p, t, ~, ~, Ve, pmid, K] = GenerateMesh(xn, yn, h, BDY, D, []) ;
% INITIAL DESIGN VARIABLES
xupp = ones(length(t), 1); xlow = ones(length(t), 1)*xmin;
x = max(xlow,min(xupp,ones(length(t),1)*vol_con));
clf ; colormap summer ; patch( 'Faces', t,'Vertices', p, 'FaceVertexCData', x, 'FaceColor', 'flat') ; colorbar ;
% PREPARE FILTER
[ H, Hs] = FilterIndex(pmid, rmin, Ve);  % PREPARE FILTER
%% MAIN LOOP
loop = 0 ; loopbeta = 0; beta = 1e-6; tol = 0;
while tol == 0
    xold = x; loop = loop + 1 ; loopbeta = loopbeta + 1;
    % FEA AND SENSITIVITY ANALYSIS
    [dCe, J0] = FEA(t, p, BDY, x, K) ;
    dVe = Ve./sum(Ve); % SENSITIVITY OF VOLUME
    vol = sum( Ve.* x)/( lelx * lely) ;
    % GRAYNESS DEGREE
    WetE = Ve .* min(abs(x - 0), abs(x - 1));
    delta = sum(WetE) / sum(Ve);
    % OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
    [x] = OcUpdate(xold, dCe, dVe, Ve, vol_con, lelx, lely, H, Hs, beta, xupp, xlow) ;
    % RESULT DISPLAY
    clf; colormap summer ;
    patch( 'Faces', t, 'Vertices', p, 'FaceVertexCData', x, 'FaceColor', 'flat') ;
    colorbar ; hold on ; axis off equal tight ; drawnow ;
    change = max(abs(xold(:)-x(:))) ;
    fprintf( 'It.:%5i Obj.:%8.4f Vol.:%7.3f ch.:%7.3f 0/1:%7.3f\n', loop, J0, vol, change, delta) ;
    % INCREASE PROJECTION PARAMETER
    if ((loopbeta >= 5 && change < 0.02) || loopbeta >= 50)
        if delta < 0.02 && change < 0.02
            tol = 1;
        end
        if delta >= 0.02
            beta = beta + dbeta; loopbeta = 0;
            fprintf('Parameter beta increased to %g.\n',beta) ;
            if delta < 0.1
                xBF = griddata(pmid(:,1), pmid(:,2), x , xn, yn, 'linear') ;
                xBD = griddata(pmid(:,1), pmid(:,2), x , xn, yn, 'nearest') ;
                xBF(isnan(xBF)) = xBD(isnan(xBF)) ;
                [c] = ContourPoints(contour(xn, yn, xBF, [0.5 0.5]), d1, d2) ;
                [p,t,t1,t2,Ve,pmid,K] = GenerateMesh(xn, yn, h, BDY, D, c, xBF, Meshdiff) ;    % REMESHING THE DESIGN DOMAIN
                patch('Faces', t1, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[0 127 102]/255) ;  axis off equal tight
                patch('Faces', t2, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[255 255 102]/255) ; set(gca,'YDir','reverse') ;
                x = ones(length(t), 1); x(1:length(t1)) = xmin ;
                xupp = ones(length(t), 1); xlow = ones(length(t), 1)*xmin ;
                [H, Hs] = FilterIndex(pmid, rmin, Ve); % PREPARE FILTER FOR NEW MESH
            end
        end
    end
end
%% FINAL RESULT PLOTTTING
xBF = griddata(pmid(:,1), pmid(:,2), x , xn, yn, 'linear') ;
xBD = griddata(pmid(:,1), pmid(:,2), x , xn, yn, 'nearest') ;
xBF(isnan(xBF)) = xBD(isnan(xBF));
[c] = ContourPoints(contour(xn, yn, xBF, [0.5 0.5]), d1, d2) ;
[p,~,t1,t2] = GenerateMesh(xn, yn, h, BDY, D, c, xBF, Meshdiff) ;
patch('Faces', t1, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[0 127 102]/255) ;  axis off equal tight
patch('Faces', t2, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[255 255 102]/255) ; set(gca,'YDir','reverse') ;
end

%% FINITE ELEMENT ANALYSIS
function [dCe, J] = FEA(t, p, BDY, x, K)
NT = length(t); NDOF = 2*length(p);
elemDof = zeros(NT,6); elemDof(:,[1 3 5]) = 2*t-1; elemDof(:,[2 4 6]) = 2*t;
iK = reshape( kron( elemDof, ones(6,1))', 36*NT,1) ;
jK = reshape( kron( elemDof, ones(1,6))', 36*NT,1) ;
xi = repelem( x, 36, 1) ;
sK = xi.*reshape( K, 36*NT, 1) ;
NK = fsparse( iK, jK, sK, [NDOF, NDOF]) ;
fixedNodes = find(p(:,1)==BDY(1,1)) ;
forceNodes = find(p(:,1)==BDY(2,1) & p(:,2)==(BDY(1,2)+BDY(2,2))/2) ;
fixedDof = [2*fixedNodes-1; 2*fixedNodes] ;
AllDof = 1:NDOF ;
freeDofs = setdiff( AllDof, fixedDof) ;
F = sparse( 2*forceNodes, 1, -1, 2*length(p), 1) ;
U = zeros(NDOF, 1) ;
GK = decomposition(NK(freeDofs,freeDofs),'chol','lower');
U(freeDofs) = GK \ F(freeDofs) ;
J = F' * U / 2;
dCe = deal(zeros(NT,1));
for i = 1 : NT
    dCe(i) = - U(elemDof(i,:))'*K(:,6*i-5:6*i)*U(elemDof(i,:)) ;
end
end

%% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
function [xnew] = OcUpdate(xold, dCe, dVe, Ve, vol_con, lelx, lely, H, Hs, beta, xupp, xlow)
l1 = 0 ; l2 = 100 ;  move = 0.2 ;
while (l2-l1)/(l1+l2) > 1.0e-4
    lmid = 0.5*(l2+l1) ;
    x = max(xlow, max( xold-move, min( xupp, min( xold+move, xold.*sqrt(max(0,-lmid.*dCe./dVe)))))) ;
    xfilt = max(xlow,min(xupp,H*x(:)./Hs));
    xnew = (tanh(beta*0.5)+tanh(beta*(xfilt-0.5)))/2/tanh(beta*0.5);
    xnew = max(xlow, max(xold-move, min( xupp, min(xold+move, xnew)))) ;
    if sum(xnew(:).*Ve(:)) < vol_con*lelx*lely, l1 = lmid; else l2 = lmid ; end
end
end

%% FIND POINTS ON THE CONTOUR
function [c] = ContourPoints(c,d1,d2)
num = 1; col = 1;
while col < size(c,2)
    idx = col+1:col+c(2,col);
    s(num).ck = c(:,idx);
    s(num).sed = norm(s(num).ck(:,1)-s(num).ck(:,end));
    s(num).isopen = s(num).sed>1e-12;
    num = num+1; col = col+c(2,col)+1;
end
c = [];
for k = 1:num-1
    ck = s(k).ck;
    if length(ck)>4
        ndist = vecnorm(diff(ck,1,2),2,1);
        ck(:,find(ndist < d1) +1) = 0.5*(ck(:,ndist < d1)+ck(:,find(ndist < d1) +1));  % REMOVE TOO-CLOSE NODES
        ck(:,ndist < d1) = [];
        if s(k).sed < d1
            ck(:, end) = [];
        end
        if ~isempty(ck)
            if ~s(k).isopen && ~isempty(ck)
                ck = [ck ck(:,1)];
            end
            ndist = vecnorm(diff(ck,1,2),2,1);
            ct = ck; ck = ct(:,1);
            for i = 1:length(ndist)
                if  ndist(i) > d2
                    ck = [ck 0.5*(ct(:,i)+ct(:,(i+1))) ct(:,i+1)]; % INCLUDE ADDITIONAL NODES
                else
                    ck = [ck ct(:,i+1)];
                end
            end
        end
        c = [c; ck(:,1:end-(~s(k).isopen))'];
    end
end
end

%% BODY FITTED MESH GENERATOR
function [p,t,t1,t2,Ve,pmid,K] = GenerateMesh(xn, yn, h, BDY, D, c, dN, Meshdiff)
if isempty(c) == 1
    [x,y] = meshgrid(BDY(1,1):h: BDY(2,1), BDY(1,2):h: BDY(2,2)) ;
    x(2:2:end,:)=x(2:2:end,:)+0.5*h; pi=[x(:),y(:)];
    pi(:,1) = min(BDY(2,1),max(BDY(1,1),pi(:,1))); pi(:,2) = min(BDY(2,2),max(BDY(1,2),pi(:,2)));
    p = [ BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2; pi];
    p = unique( p, 'rows', 'stable');
    t = delaunay(p);
    t1 = [] ; t2 = [] ;
else
    [x,y] = meshgrid(BDY(1,1): h: BDY(2,1), BDY(1,2): sqrt(3)/2*h : BDY(2,2)) ;
    x(2:2:end,:)=x(2:2:end,:)+0.5*h; pi=[x(:),y(:)];
    % d = min(max(2*Emin,min(pdist2(pi, c), [], 2)),Meshdiff);
    dist = min(pdist2(pi, c), [], 2);
    d = rescale(dist, h, Meshdiff*h);
    r0 = 1./d.^2;
    pfix=[c; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2];
    pfix = unique(pfix,'rows','stable');
    load('randommatrixmc.mat','ranmatrix');   % LOAD SAVED RANDOM MATRIX
    % ranmatrix = rand(size(pi,1),1);       % OR GENERATE A NEW RANDOM MATRIX
    % save('randommatrixmc.mat','ranmatrix');   % SAVE RANDOM MATRIX
    p=[pfix; pi(ranmatrix(1:size(pi,1))<r0./max(r0),:)];  % DELETE NODES USING REJECTION METHOD
    p(:,1) = min(BDY(2,1),max(BDY(1,1),p(:,1))); p(:,2) = min(BDY(2,2),max(BDY(1,2),p(:,2)));
    L1xy = reshape(sqrt(min((xn(:) - c(:,1)').^2 + (yn(:) - c(:,2)').^2, [], 2)),size(xn));
    L1xy = rescale(L1xy, h, Meshdiff*h); % DESIRED EDGE LENGTH L1
    p1 = 0;
    % NOVE-MOVING LOOPS
    for i = 1:200
        if max(sum((p-p1).^2,2))>0.1*h
            p = unique(p,'rows','stable');
            t = delaunay(p);  % DELAUNAY TRIANGULATION
            edges = unique(sort([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],2),'rows') ;
            p1 = p;
        end
        L = sqrt(sum((p(edges(:,1),:)-p(edges(:,2),:)).^2,2));
        pm = (p(edges(:,1),:)+p(edges(:,2),:))/2;
        L1 = interp2(xn,yn,L1xy,pm(:,1),pm(:,2),'cubic');  % INTERPOLATE DENSITY INTO EDGE CENTROIDS
        L0 = 1.2*L1*sqrt(sum(L.^2)/sum(L1.^2));
        Fb = max(L0-L,0)./L *[1,1].*(p(edges(:,1),:)-p(edges(:,2),:));
        Fp = full(sparse(edges(:,[1,1,2,2]),ones(size(L1))*[1,2,1,2],[Fb,-Fb],size(p,1),2));
        Fp(1:size(pfix,1),:) = 0 ;
        p = p + 0.2*Fp ;  % MOVE NODE ACCORDING TO VIRTUAL FORCE
        p(:,1) = min(BDY(2,1),max(BDY(1,1),p(:,1))) ;
        p(:,2) = min(BDY(2,2),max(BDY(1,2),p(:,2))) ;
        if max(sum(0.2*Fp,2))<0.05*h, break; end
    end
    Ve = zeros(length(t),1) ;
    for i = 1:length(t)
        J = [p(t(i,1),1)-p(t(i,3),1) p(t(i,1),2)-p(t(i,3),2) ; p(t(i,2),1)-p(t(i,3),1) p(t(i,2),2)-p(t(i,3),2)];
        Ve(i) = 0.5.*det(J);  % ELEMENT VOLUME
    end
    t(Ve == 0,:) = [];
    pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
    dE = interp2(xn,yn,dN,pmid(:,1),pmid(:,2),'cubic');  % INTERPOLATE DENSITY INTO ELEMENT CENTROIDS
    t1=[t(dE<0.5,:)]; t2=t(dE>=0.5,:); t=[t1;t2];
end
pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
Ve = zeros(length(t),1) ;
for i = 1:length(t)
    J = [p(t(i,1),1)-p(t(i,3),1) p(t(i,1),2)-p(t(i,3),2) ; p(t(i,2),1)-p(t(i,3),1) p(t(i,2),2)-p(t(i,3),2)];
    Ve(i) = 0.5.*det(J);  % ELEMENT VOLUME
    Be = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0; 0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1); -J(2,1) J(2,2) J(1,1) -J(1,2) J(2,1)-J(1,1) -J(2,2)+J(1,2)];
    K(:,6*i-5:6*i) = det(J)/2*Be'*D*Be;
end
end

%% PREPARE FILTER
function [H, Hs] = FilterIndex(pmid,rmin,Ve)
maxN =  ceil (100 * rmin);    NT = length (pmid);
iH = ones(NT*maxN, 1) ;
jH = ones(NT*maxN, 1) ;
sH = zeros(NT*maxN, 1) ;
Start = 1;
for i = 1 : NT
    pni = find((pmid(:,1) - pmid(i,1)).^2 + (pmid(:,2) - pmid(i,2)).^2 <= rmin^2);
    Weight = max((rmin - pdist2(pmid(i,:),pmid(pni,:))),0) .* (Ve(pni))';
    pEnd = Start + length(Weight);
    iH(Start:pEnd-1) = i;
    jH(Start:pEnd-1) = pni;
    sH(Start:pEnd-1) = Weight;
    Start = pEnd;
end
iH = iH(1:Start - 1);
jH = jH(1:Start - 1);
sH = sH(1:Start - 1);
H = sparse( iH, jH, sH) ;
Hs = sum( H, 2) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by ZC. Zhuang and Y. M. Xie                %
% Centre for Innovative Structures and Materials, RMIT University         %
% Please send your comments to: zhuanginhongkong@outlook.com              %
%                                                                         %
% The program is proposed for educational purposes and is introduced in   %
% the paper - An Integrated Framework for Body-Fitted Topology            %
% Optimization Considering Buckling Load Factor, Stiffness,               %
% and von Mises Stress, CMAME, 2025                                       %
%                                                                         %
% Disclaimer:                                                             %
% The authors reserve all rights but do not guarantee that the code is    %
% free from errors. Furthermore, we shall not be liable in any event.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%