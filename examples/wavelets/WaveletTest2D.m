N = 256;
impath = '../../StandardTestImages/';
ln = imread([impath 'Lena512.png'],'png'); ln = double(imresize(ln,[N N]));
hs = imread([impath 'house.png'],'png'); hs = double(imresize(hs,[N N]));
pt = imread([impath 'Cameraman256.png'],'png'); pt = double(imresize(pt,[N N]));
%ot = imread([impath 'ot.pgm']); ot = double(imresize(ot,[N,N]));
HG = fspecial('gaussian',N,N/32);
sl = phantom(N);
Deblur = 0; % set to 0 for denoising
xt =  pt; %pt;%hs;%sl
xt = xt./max(max(xt)); % rescale to 1;
Tf = 0.1; % fraction to apply to max to get threshold
figure(1); clf;
subplot(2,3,1); colormap(gray); imagesc(xt); title('original');
immax = max(max(xt));
%% ------------------------- blur image ----------------------------------

sig = 0.02;
noiselevel = 0.1;
if Deblur 
xb = imblur(xt,sig);
else
    sig = 0;
    xb = xt;
end
xb = xb + noiselevel*immax.*randn(size(xt));
figure(1); 
subplot(2,3,4); colormap(gray); imagesc(xb); title(['blurred +',num2str(100*noiselevel),'% noise']);
%% --------------- wavelet decomposition ---------------------------------
NW = 3;
%wfilt = 'haar'; 
wfilt = 'db4';
[C,S] = wavedec2(xt,NW,wfilt);
CA = zeros(S(1,:));
cind = 1;
cend = length(reshape(CA,[],1));
CA = reshape(C(cind:cend),S(1,:));
Cim = zeros(N,N);
Cim(1:S(1,1),1:S(1,2)) = CA./max(C);
for n = 1:NW
    sn = S(n,:);
    CH{n} = zeros(S(n+1,:));
    len = length(reshape(CH{n},[],1));
    disp(['level ',num2str(n), ' coeff length ',num2str(len)]);
    cind = cend + 1;
    cend = cind + len-1;
    CH{n} = reshape(C(cind:cend),S(n+1,:));
    Cim(S(n+1,1)+1:2*S(n+1,1),1:S(n+1,:)) = CH{n}/max(C);
    cind = cend + 1;
    cend = cind + len-1;
    CV{n} = reshape(C(cind:cend),S(n+1,:));
    Cim(1:S(n+1,1),S(n+1,1)+1:2*S(n+1,:)) = CV{n}/max(C);    
    cind = cend + 1;
    cend = cind + len-1;
    CD{n} = reshape(C(cind:cend),S(n+1,:));
    Cim(S(n+1,1)+1:2*S(n+1,1),S(n+1,1)+1:2*S(n+1,:)) = CD{n}/max(C);        
end
figure(2); imagesc(Cim); colormap(gray); colorbar;

%% look at wavelets
for k = 64:100:128*128
CT = zeros(size(C));
CT(k) = 1;
xtt = waverec2(CT,S,wfilt);
figure(5); clf; imagesc(xtt);
pause(0.25);
end

%% - test setting wavelet coefficients to zero
CR = C;
cz = find(abs(C)<0.01*max(C));
CR(cz) = 0;
cfact = length(cz)/length(C);
disp(['setting ',num2str(100*cfact),'% of coefficient to zero']);
xr = waverec2(CR,S,wfilt);
figure(1);
subplot(2,3,2); colormap(gray); imagesc(xr); title(['original, ',wfilt,' wavelet compressed ',num2str(100*cfact),'%']);
subplot(2,3,3); colormap(gray); imagesc(xr-xt); title(['difference']);

%% -------- wavelet decomposition - noisy case ----------------------------
[Cb,Sb] = wavedec2(xb,NW,wfilt);
CA = zeros(S(1,:));
cind = 1;
cend = length(reshape(CA,[],1));
CA = reshape(C(cind:cend),S(1,:));
Cim = zeros(N,N);
Cim(1:S(1,1),1:S(1,2)) = CA./max(C);
for n = 1:NW
    sn = S(n,:);
    CH{n} = zeros(S(n+1,:));
    len = length(reshape(CH{n},[],1));
    disp(['level ',num2str(n), ' coeff length ',num2str(len)]);
    cind = cend + 1;
    cend = cind + len-1;
    CH{n} = reshape(C(cind:cend),S(n+1,:));
    Cim(S(n+1,1)+1:2*S(n+1,1),1:S(n+1,:)) = CH{n}/max(C);
    cind = cend + 1;
    cend = cind + len-1;
    CV{n} = reshape(C(cind:cend),S(n+1,:));
    Cim(1:S(n+1,1),S(n+1,1)+1:2*S(n+1,:)) = CV{n}/max(C);    
    cind = cend + 1;
    cend = cind + len-1;
    CD{n} = reshape(C(cind:cend),S(n+1,:));
    Cim(S(n+1,1)+1:2*S(n+1,1),S(n+1,1)+1:2*S(n+1,:)) = CD{n}/max(C);        
end
figure(3);clf;
imagesc(Cim); colormap(gray); colorbar;
figure(4);clf;
imagesc(log(abs(Cim))); colormap(gray); colorbar;

%% - test setting wavelet coefficients to zero
CRb = Cb;
cz = find(abs(Cb)<0.025*max(Cb));
CRb(cz) = 0;
cfact = length(cz)/length(Cb);
disp(['setting ',num2str(100*cfact),'% of coefficient to zero']);
xrb = waverec2(CRb,Sb,wfilt);
figure(1);
subplot(2,3,5); colormap(gray); imagesc(xrb); title(['blurred, ',wfilt,' wavelet compressed ',num2str(100*cfact),'%']);
subplot(2,3,6); colormap(gray); imagesc(xrb-xb); title('difference');
