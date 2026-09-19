function [x,Zr,Zc,nOut] = decodeMixedDesign(y,p)
%DECODEMIXEDDESIGN Decode continuous variables and integer set indices.

validateattributes(y,{'numeric'},{'real','vector','finite'});
y = y(:).';
if numel(y)~=10
    error('CycloidalDrive:InvalidMixedVector','Mixed design vector must have 10 values.');
end
iZr = round(y(9));
iOut = round(y(10));
if iZr<1 || iZr>numel(p.Zr_set) || iOut<1 || iOut>numel(p.nOut_set)
    error('CycloidalDrive:DiscreteIndexOutOfRange','Discrete design index is out of range.');
end
x = y(1:8);
Zr = p.Zr_set(iZr);
Zc = Zr-1;
nOut = p.nOut_set(iOut);
end
