function f = objectiveMixedWrapper(y,p)
[x,Zr,Zc,nOut] = decodeMixedDesign(y,p);
f = objectiveWrapper(x,Zr,Zc,nOut,p);
end
