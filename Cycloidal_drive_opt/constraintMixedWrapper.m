function [c,ceq] = constraintMixedWrapper(y,p)
[x,Zr,Zc,nOut] = decodeMixedDesign(y,p);
[c,ceq] = constraintWrapper(x,Zr,Zc,nOut,p);
end
