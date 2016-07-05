%% Bipul Mohanto
% Color in Informatics and MEdia Technology
% contact: bipul.mohanto@yahoo.com
% Algorithm Implimentation: 1. Template Matching, 2. Texture, 3. Quadtree Decomposition, 4. Descriptors

clear all
close all
clc
 
%% template matching
I=imread('p1.bmp');
I=rgb2gray(I);
 
Ob=imread('letter.bmp');
Ob = rgb2gray(Ob);
 
cc = normxcorr2(Ob,I);
% best score where best match
[max_cc, imax] = max(abs(cc(:)));
[ypeak, xpeak] = ind2sub(size(cc),imax(1));
point = [ (ypeak-size(Ob,1)/2+1) (xpeak-size(Ob,2)/2+1) ];


K=cat(3,I,I,I);
if (max_cc>0.8)
    K(point(1), point(2),2)=255; % good correlation in green
else
    K(point(1), point(2),1)=255; % not so good correlation in red
end
figure, imshow(K);

% results : match 100% for p1, 78% for p2 and 84% for p3
% 80% here seems to be a good threshlod to select only letters T

%% co-occurrence matrix
clc
clear all
close all

I=imread('fur_txtr.tif');
I2=imread('fur_txtr_noisy.tif');
for d=1:50 % horizontal direction, 
    glcm = graycomatrix(I,'NumLevels',256,'Offset',[0 d],'Symmetric',true);
    stats = graycoprops(glcm,{'contrast','homogeneity','energy','correlation'});
    tabICont(d)=stats.Contrast;
    tabIHomo(d)=stats.Homogeneity;
    tabIEner(d)=stats.Energy;
    tabICorr(d)=stats.Correlation;
    
    glcm = graycomatrix(I2,'NumLevels',256,'Offset',[0 d],'Symmetric',true);
    stats = graycoprops(glcm,{'contrast','homogeneity','energy','correlation'});
    tabI2Cont(d)=stats.Contrast;
    tabI2Homo(d)=stats.Homogeneity;
    tabI2Ener(d)=stats.Energy;
    tabI2Corr(d)=stats.Correlation;

end

ecart = abs(tabICorr-tabI2Corr);
plot(ecart); % to show the good robustness against noise notably when the d distance increases

%% quadtree decompostion
I = imread('cheetah.bmp');

S = qtdecomp(I,0.5);
blocks = repmat(uint8(0),size(S));
for dim = [256 128 64 32 16 8 4 2 1];    
  numblocks = length(find(S==dim));    
  if (numblocks > 0)        
    values = repmat(uint8(1),[dim dim numblocks]);
    values(2:dim,2:dim,:) = 0;
    blocks = qtsetblk(blocks,S,dim,values);
  end
end

blocks(end,1:end) = 1;
blocks(1:end,end) = 1;

imshow(I), figure, imshow(blocks,[])

%% cards and signatures

clear all
close all
clc

I=imread('cards.bmp');

BW = im2bw(I);
figure(1), imshow(BW);

L = bwlabel(BW);
s  = regionprops(L,'Solidity');

% heart
idx = find([s.Solidity] > 0.92);
BW2 = ismember(L,idx);
outCo = BW2.*1;

% diamond
idx = find([s.Solidity] > 0.88 & [s.Solidity] <= 0.92);
BW2 = ismember(L,idx);
outCa = BW2.*2;

% spade
idx = find([s.Solidity] > 0.80 & [s.Solidity] <= 0.88);
BW2 = ismember(L,idx);
outPi = BW2.*3;

% club
idx = find([s.Solidity] <= 0.80);
BW2 = ismember(L,idx);
outTr = BW2.*4;

L = outCo+outCa+outPi+outTr;
RGB = label2rgb(L,'jet',[0 0 0]); 
figure, imshow(RGB);

