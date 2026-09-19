function [c,ceq] = constraintWrapper(x,Zr,Zc,nOut,p)
%CONSTRAINTWRAPPER Dimensionless, order-one nonlinear inequalities.

d = evaluateDesignCached(x,Zr,Zc,nOut,p);
if ~d.isValid
    % Individual margins are retained where finite; the first entry provides
    % an unambiguous order-one invalid-state flag.
    invalid = 1;
else
    invalid = -1;
end

c = [ ...
    invalid;
    d.profileIntersect;
    d.outputHoleViolation/max(p.minLigamentReq,realmin);
    (p.minLigamentReq-d.minLigament)/p.minLigamentReq;
    -d.ringPinSpacingMargin/max(p.ringPinAssemblyMargin,realmin);
    -d.curvatureMargin/max(p.minCurvatureRadius,realmin);
    -d.symmetryMargin/max(1e-3*d.Rp,realmin);
    -d.shortWidthMargin/max(p.shortWidthCoefficientMin,realmin);
    (d.maxRingGap-p.maxNoLoadRingGap)/p.maxNoLoadRingGap;
    (-p.maxRingInterference-d.minRingGap)/p.maxRingInterference;
    (p.tDisc_min_req-d.tDisc)/p.tDisc_min_req;
    (p.SFmin_bend-d.SF_bend)/p.SFmin_bend;
    (p.SFmin_bearing-d.SF_bearing)/p.SFmin_bearing;
    (p.SFmin_shear-d.SF_shear)/p.SFmin_shear;
    (p.SFmin_pin-d.SF_pin)/p.SFmin_pin;
    (p.SFmin_creep-d.SF_creep)/p.SFmin_creep;
    (p.SFmin_wear-d.SF_wear)/p.SFmin_wear;
    (d.backlash_deg-p.backlashMax_deg)/p.backlashMax_deg;
    d.torqueBalanceError_Nm/max(abs(p.Tout_nominal),realmin)-1e-6
    ];
c(~isfinite(c)) = 1;
ceq = [];
end
