function d = evaluateDesign(x,Zr,Zc,nOut,p)
%EVALUATEDESIGN Pure deterministic nominal design evaluation.

v = unpackDesign(x,p);
g = cycloidGeometry(v.Rp,v.e,v.dRing,v.dOut,v.tDisc, ...
    v.ringClearance,v.rOutPitch,Zr,Zc,nOut,p,v.outputClearance);
a = evaluatePhaseSweep(g,p,Zc);
s = a.s;
b = a.b;

d.Zr = Zr;
d.Zc = Zc;
d.reductionRatio = Zc;
d.nOut = nOut;
d.x = v.x;
d.Rp = v.Rp;
d.e = v.e;
d.dRing = v.dRing;
d.dOut = v.dOut;
d.tDisc = v.tDisc;
d.clearance = v.ringClearance; % Legacy alias.
d.ringClearance = v.ringClearance;
d.outputClearance = v.outputClearance;
d.rHole = v.rOutPitch; % Legacy alias.
d.rOutPitch = v.rOutPitch;

d.mass_kg = g.mass_kg;
d.backlash_deg = b.backlash_deg;
d.loadedLostMotion_deg = b.loaded_lost_motion_deg;
d.kt_Nm_per_rad = b.kt_Nm_per_rad;
d.sigmaContact_MPa = s.sigmaContact_MPa;
d.sigmaBearingOut_MPa = s.sigmaBearingOut_MPa;
d.sigmaRoot_MPa = s.sigmaRoot_MPa;
d.tauWeb_MPa = s.tauWeb_MPa;
d.sigmaPinRing_MPa = s.sigmaPinRing_MPa;
d.sigmaPinOut_MPa = s.sigmaPinOut_MPa;

d.SF_bend = s.SF_bend;
d.SF_bearing = s.SF_bearing;
d.SF_shear = s.SF_shear;
d.SF_pin = s.SF_pin;
d.SF_creep = s.SF_creep;
d.SF_wear = s.SF_wear;

d.profileIntersect = g.profileIntersect;
d.outputHoleViolation = g.outputHoleViolation;
d.minLigament = g.minLigament;
d.ringPinSpacingMargin = g.ringPinSpacingMargin;
d.curvatureMargin = g.curvatureMargin;
d.symmetryMargin = g.symmetryMargin;
d.minRingGap = a.minMinRingGap;
d.maxRingGap = a.maxMinRingGap;
d.shortWidthCoefficient = g.shortWidthCoefficient;
d.shortWidthMargin = g.shortWidthMargin;
d.contactCountRing = a.minRingContactCount;
d.contactCountOutput = a.minOutputContactCount;
d.bearingReaction_N = a.maxBearingReaction_N;
d.torqueBalanceError_Nm = a.maxTorqueBalanceError;
d.evaluatedInputPhases_rad = a.phases;

d.marginMin = min([ ...
    (d.SF_bend-p.SFmin_bend)/p.SFmin_bend, ...
    (d.SF_bearing-p.SFmin_bearing)/p.SFmin_bearing, ...
    (d.SF_shear-p.SFmin_shear)/p.SFmin_shear, ...
    (d.SF_pin-p.SFmin_pin)/p.SFmin_pin, ...
    (d.SF_creep-p.SFmin_creep)/p.SFmin_creep, ...
    (d.SF_wear-p.SFmin_wear)/p.SFmin_wear]);

d.isValidGeometry = g.isValidGeometryBase && a.allGeometryValid;
d.isValidContact = a.allContactsValid;
d.isValidStress = s.isValidStress;
d.isValidBacklash = b.isValidBacklash;
d.isValid = d.isValidGeometry && d.isValidContact && ...
    d.isValidStress && d.isValidBacklash && isfinite(d.mass_kg);

% Robust values are populated only in the post-optimization validation pass.
d.robustFailureProb = NaN;
d.robustFailureProbLower95 = NaN;
d.robustFailureProbUpper95 = NaN;
d.robustMeanBacklash_deg = NaN;
d.robustMeanMinSF = NaN;
end
