function [Reference, CorrectedRef,tform] = maskShift(Reference,ImageNew)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

ImageNew = double(ImageNew);
[Corrected, CorrectedRef,tform] = registerImages(Reference,ImageNew);
Corrected.Image = ImageNew;
Reference = Corrected;
end

