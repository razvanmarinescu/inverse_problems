function compare_regselect(K,RH,f,yn,sigma,figno);
%
% test different regularisation parameter selection methods.
%
nsol = size(f,1); % size of solution
%% ------------- Discrepancy Principle --------------
logalpha = [-3:0.1:0];
alpha = 10.^logalpha;
%alpha = [0:0.001:0.1];
dp = zeros(size(alpha));
ed = dp;
for k = 1:size(alpha,2)
    [dp(k),edp(k)] = DP(yn,K,f,sigma,alpha(k),RH);
end;
figure(5);
semilogx(alpha,dp); title('Discrepancy Principle value');
%hold on; semilogx(alpha,log(edp),'ko');

% binary search in DP
lo = alpha(1); hi = alpha(end); 
dplo = DP(yn,K,f,sigma,lo,RH);
dphi = DP(yn,K,f,sigma,hi,RH);
dpmid = dphi;
dptol = 1e-6*dphi;
while(abs(dpmid) > dptol)
    mid = (lo+hi)/2;
    dpmid = DP(yn,K,f,sigma,mid,RH);
    if(dpmid < 0)
        lo = mid;
    else
        hi = mid;
    end
end
lambdaDP = mid;
[dpopt,edpopt] = DP(yn,K,f,sigma,lambdaDP,RH);
disp(['DP - \alpha=',num2str(lambdaDP), ' DP = ',num2str(dpopt)]);
% display this
fDP = (K'*K + lambdaDP*RH)\(K'*yn);
figure(figno); clf; hold on;
plot(f,'r');title('comparison of regularisation selection methods');
plot(fDP,'--m+');
legend('True solution',['DP - \alpha=',num2str(lambdaDP)]);

%% ------------- Miller Criterion --------------
logalpha = [-3:0.1:0];
alpha = 10.^logalpha;
%alpha = [0:0.001:0.1];
miller = zeros(size(alpha));
em = miller;
for k = 1:size(alpha,2)
    [miller(k),em(k)] = Miller(yn,K,f,sigma,alpha(k),RH);
end;
figure(6);
semilogx(alpha,miller); title('Miller Criterion');
hold on; semilogx(alpha,log(em),'ko');

% binary search in Miller
lo = alpha(1); hi = alpha(end); 
mlo = Miller(yn,K,f,sigma,lo,RH);
mhi = Miller(yn,K,f,sigma,hi,RH);
mmid = mhi;
mtol = 1e-6*mhi;
while(abs(mmid) > mtol)
    mid = (lo+hi)/2;
    mmid = Miller(yn,K,f,sigma,mid,RH);
    if(mmid < 0)
        lo = mid;
    else
        hi = mid;
    end
end
lambdaMiller = mid;
disp(['Miller - \alpha=',num2str(lambdaMiller), ' Miller = ',num2str(mmid)]);
% display this
fm = (K'*K + lambdaMiller*RH)\(K'*yn);
figure(figno); 
plot(fm,'--c+');
legend('True solution',['DP - \alpha=',num2str(lambdaDP)],['Miller - \alpha=',num2str(lambdaMiller)]);

%% ------------- Predictive Risk --------------
logalpha = [-3:0.1:0];
alpha = 10.^logalpha;
%alpha = [0:0.001:0.1];
p = zeros(size(alpha));
ep = p;
for k = 1:size(alpha,2)
    [p(k),ep(k)] = predrisk(yn,K,f,sigma,alpha(k),RH);
end;
figure(7);
loglog(alpha,p); title('Predictive Risk value');
hold on; loglog(alpha,ep,'ko');

% binary search in predrisk (this is search for minumum)
lo = alpha(1); hi = alpha(end); 
plo = predrisk(yn,K,f,sigma,lo,RH);
phi = predrisk(yn,K,f,sigma,hi,RH);
pmid = phi;
ptol = 1e-6*phi;
while(abs(pmid-plo) > ptol)
    mid = (lo+hi)/2;
    pmid = predrisk(yn,K,f,sigma,mid,RH);
    a = [lo,mid,hi];
    v = [plo,pmid,phi];
    [pv,perm]=sort(v);
    lo = a(perm(1));plo = v(perm(1));
    hi = a(perm(2));phi = v(perm(2));
    pmid = v(perm(3));
end
lambdapred = mid;
disp(['PredRisk - \alpha=',num2str(lambdapred), ' PredRisk = ',num2str(pmid)]);
figure(7);loglog(lambdapred,pmid,'r+');
% display this
fpred = (K'*K + lambdapred*RH)\(K'*yn);
figure(figno); 
plot(fpred,'--k+');
legend('True solution',['DP - \alpha=',num2str(lambdaDP)],['Miller - \alpha=',num2str(lambdaMiller)],['pred - \alpha=',num2str(lambdapred)]);


%% ------------- Unbiased Predictive Risk Estimator --------------
logalpha = [-3:0.1:0];
alpha = 10.^logalpha;
%alpha = [0:0.001:0.1];
upre = zeros(size(alpha));
eu = upre;
for k = 1:size(alpha,2);
    [upre(k),eu(k)] = UPRE(yn,K,f,sigma,alpha(k),RH);
end;
figure(8);
loglog(alpha,upre); title('Unbiased Predictive Risk Estimator');
hold on; loglog(alpha,eu,'ko');

% binary search in UPRE (this is search for minumum)
lo = alpha(1); hi = alpha(end); 
ulo = UPRE(yn,K,f,sigma,lo,RH);
uhi = UPRE(yn,K,f,sigma,hi,RH);
umid = uhi;
utol = 1e-6*uhi;
while(abs(umid-ulo) > utol)
    mid = (lo+hi)/2;
    umid = UPRE(yn,K,f,sigma,mid,RH);
    a = [lo,mid,hi];
    v = [ulo,umid,uhi];
    [pv,perm]=sort(v);
    lo = a(perm(1));ulo = v(perm(1));
    hi = a(perm(2));uhi = v(perm(2));
    umid = v(perm(3));
end
lambdaUPRE = mid;
disp(['UPRE - \alpha=',num2str(lambdaUPRE), ' UPRE = ',num2str(umid)]);
figure(8);loglog(lambdaUPRE,umid,'r+');
% display this
fupre = (K'*K + lambdaUPRE*RH)\(K'*yn);
figure(figno); 
plot(fupre,'--g+');
legend('True solution',['DP - \alpha=',num2str(lambdaDP)],['Miller - \alpha=',num2str(lambdaMiller)],['pred - \alpha=',num2str(lambdapred)],['UPRE - \alpha=',num2str(lambdaUPRE)]);

%% ------------- Generalised Cross -Validation Index --------------
logalpha = [-3:0.1:0];
alpha = 10.^logalpha;
%alpha = [0:0.001:0.1];
gcv = zeros(size(alpha));
eg = gcv;
nf = zeros(size(alpha));
for k = 1:size(alpha,2);
    [gcv(k),eg(k),nf(k)] = GCV(yn,K,f,sigma,alpha(k),RH);
end;
figure(9);
loglog(alpha,gcv); title('Generalised Cross-Validation');
hold on; loglog(alpha,eg,'ko');

% binary search in GCV (this is search for minumum)
lo = alpha(1); hi = alpha(end); 
glo = GCV(yn,K,f,sigma,lo,RH);
ghi = GCV(yn,K,f,sigma,hi,RH);
gmid = ghi;
gtol = 1e-6*ghi;
while(abs(gmid-glo) > gtol)
    mid = (lo+hi)/2;
    gmid = GCV(yn,K,f,sigma,mid,RH);
    a = [lo,mid,hi];
    v = [glo,gmid,ghi];
    [pv,perm]=sort(v);
    lo = a(perm(1));glo = v(perm(1));
    hi = a(perm(2));ghi = v(perm(2));
    gmid = v(perm(3));
end
lambdaGCV = mid;
[gopt,egopt] = GCV(yn,K,f,sigma,lambdaGCV,RH);
disp(['GCV - \alpha=',num2str(lambdaGCV), ' GCV = ',num2str(gopt)]);
figure(9);loglog(lambdaGCV,gmid,'r+');
% display this
fgcv = (K'*K + lambdaGCV*RH)\(K'*yn);
figure(figno); 
plot(fgcv,'--b+');
legend('True solution',['DP - \alpha=',num2str(lambdaDP)],['Miller - \alpha=',num2str(lambdaMiller)],['pred - \alpha=',num2str(lambdapred)],['UPRE - \alpha=',num2str(lambdaUPRE)],['GCV - \alpha=',num2str(lambdaGCV)]);

%----------------- compare --------------------------------------
figure;
loglog(alpha,upre,'k'); 
hold on;
loglog(alpha,p,'--b');
loglog(alpha,gcv,'-.m');
loglog(alpha,eg./nsol,'ko');
%% ----------------- L-curve --------------------------------------
figure;
ndat = length(yn);
%loglog(edp,ndat*(dp+sigma^2),'--r+'); title('L-curve');
logalpha = [-3:0.1:0];
alpha = 10.^logalpha;
for k = 1:size(alpha,2)
    fk = (K'*K + alpha(k)*RH)\(K'*yn);
    llhd(k) = 0.5*norm(yn - K*fk)^2;
    pr(k) = 0.5*(fk)'*RH*(fk);
end;
loglog(pr,llhd,'r'); title('L-curve');
%% show other results on top of L-curve
hold on; 
% loglog(edpopt,dpopt+sigma^2,'m+');
% [gopt,egopt] = DP(yn,K,f,sigma,lambdaGCV,RH);
% loglog(egopt,gopt+sigma^2,'b+');
% [mopt,emopt] = DP(yn,K,f,sigma,lambdaMiller,RH);
% loglog(emopt,mopt+sigma^2,'c+');
% [popt,epopt] = DP(yn,K,f,sigma,lambdapred,RH);
% loglog(epopt,popt+sigma^2,'k+');
% [uopt,euopt] = DP(yn,K,f,sigma,lambdaUPRE,RH);
% loglog(euopt,uopt+sigma^2,'g+');

fk = (K'*K + lambdaDP*RH)\(K'*yn);
de = 0.5*norm(yn - K*fk)^2;
pe = 0.5*(fk)'*RH*(fk);
loglog(pe,de,'m+');

fk = (K'*K + lambdaGCV*RH)\(K'*yn);
de = 0.5*norm(yn - K*fk)^2;
pe = 0.5*(fk)'*RH*(fk);
loglog(pe,de,'b+');

fk = (K'*K + lambdaMiller*RH)\(K'*yn);
de = 0.5*norm(yn - K*fk)^2;
pe = 0.5*(fk)'*RH*(fk);
loglog(pe,de,'c+');

fk = (K'*K + lambdapred*RH)\(K'*yn);
de = 0.5*norm(yn - K*fk)^2;
pe = 0.5*(fk)'*RH*(fk);
loglog(pe,de,'k+');

fk = (K'*K + lambdaUPRE*RH)\(K'*yn);
de = 0.5*norm(yn - K*fk)^2;
pe = 0.5*(fk)'*RH*(fk);
loglog(pe,de,'g+');
