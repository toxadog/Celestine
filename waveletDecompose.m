function W = waveletDecompose(Image,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[n,m]=size(Image);
% parameters for wavelet analysis
if size(varargin,2)==2
sc1=varargin{1};
sc2=varargin{2};
else
sc1=20;
sc2=30;
end
dsc=1;
scale = [sc1:dsc:sc2]; 
%% 1. Wavelet analysis
% Add white boundaries to avoid boundary effect caused by wavelet transform
borderSize = 101;
% ImageBorder=zeros(n+borderSize,m+borderSize);
ImageBorder=ones(n+borderSize,m+borderSize)*mean(mean(Image));
ImageBorder(1+(borderSize-1)/2:size(ImageBorder,1)-(borderSize-1)/2-1,1+(borderSize-1)/2:size(ImageBorder,2)-(borderSize-1)/2-1)=Image;
% Wavelet decomposition
cwtmexh = cwtft2(ImageBorder,'wavelet','mexh','scales',scale);
% Compute weight coeffitiens for image reconstruction (to be
% verified)
Wcoef=ones(size(cwtmexh.cfs));
for i=1:size(Wcoef,4)
Wcoef(:,:,1,i)=Wcoef(:,:,1,i)*cwtmexh.wav_norm(i)*cwtmexh.scales(i)^2/(0.25*dsc);
end
% Reconstrution of the filtered image by summation of the wavelet transform
% components with corresponding weight coeffitiens
W=-sum(cwtmexh.cfs./Wcoef,4);
% Remove white boundaries
W=W(1+(borderSize-1)/2:size(ImageBorder,1)-(borderSize-1)/2-1,1+(borderSize-1)/2:size(ImageBorder,2)-(borderSize-1)/2-1); 
end

