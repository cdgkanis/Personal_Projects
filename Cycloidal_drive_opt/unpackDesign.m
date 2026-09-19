function v = unpackDesign(x,p)
%UNPACKDESIGN Validate and name the continuous design variables.
% Seven-value legacy inputs use p.outputClear_nom.

validateattributes(x,{'numeric'},{'real','vector','finite'},mfilename,'x');
x = x(:).';
if numel(x)==7
    x = [x, p.outputClear_nom];
elseif numel(x)~=8
    error('CycloidalDrive:InvalidDesignVector', ...
        'x must have 7 legacy values or 8 current values.');
end
if any(x<=0)
    error('CycloidalDrive:NonpositiveDesign', ...
        'All design variables must be positive.');
end

v.Rp = x(1);
v.e = x(2);
v.dRing = x(3);
v.dOut = x(4);
v.tDisc = x(5);
v.ringClearance = x(6);
v.rOutPitch = x(7);
v.outputClearance = x(8);
v.x = x;
end
