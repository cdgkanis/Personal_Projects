function d = evaluateDesignCached(x,Zr,Zc,nOut,p)
%EVALUATEDESIGNCACHED Share deterministic work across objective/constraints.

persistent lastX lastZr lastZc lastNOut lastP lastDesign
hit = ~isempty(lastX) && isequal(x,lastX) && isequal(Zr,lastZr) && ...
    isequal(Zc,lastZc) && isequal(nOut,lastNOut) && isequaln(p,lastP);
if hit
    d = lastDesign;
    return
end
d = evaluateDesign(x,Zr,Zc,nOut,p);
lastX = x;
lastZr = Zr;
lastZc = Zc;
lastNOut = nOut;
lastP = p;
lastDesign = d;
end
