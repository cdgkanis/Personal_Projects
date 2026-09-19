function a = evaluatePhaseSweep(g,p,Zc)
%EVALUATEPHASESWEEP Evaluate worst cases across a full input revolution.

nPhase = max(1,round(p.nPhaseSamples));
phases = p.inputPhase+(0:nPhase-1)'*(2*pi/nPhase);
cases = repmat(struct('g',[],'ls',[],'s',[],'b',[]),nPhase,1);
for i = 1:nPhase
    gp = applyRingContactPhase(g,p,phases(i),Zc);
    ls = loadSharingModel(gp,p);
    s = stressModel(gp,ls,p);
    b = backlashModel(gp,ls,p);
    cases(i).g = gp;
    cases(i).ls = ls;
    cases(i).s = s;
    cases(i).b = b;
end

s0 = cases(1).s;
b0 = cases(1).b;
s0.sigmaContact_MPa = max(arrayfun(@(q)q.s.sigmaContact_MPa,cases));
s0.sigmaBearingOut_MPa = max(arrayfun(@(q)q.s.sigmaBearingOut_MPa,cases));
s0.sigmaRoot_MPa = max(arrayfun(@(q)q.s.sigmaRoot_MPa,cases));
s0.tauWeb_MPa = max(arrayfun(@(q)q.s.tauWeb_MPa,cases));
s0.sigmaPinRing_MPa = max(arrayfun(@(q)q.s.sigmaPinRing_MPa,cases));
s0.sigmaPinOut_MPa = max(arrayfun(@(q)q.s.sigmaPinOut_MPa,cases));
s0.SF_bend = min(arrayfun(@(q)q.s.SF_bend,cases));
s0.SF_bearing = min(arrayfun(@(q)q.s.SF_bearing,cases));
s0.SF_shear = min(arrayfun(@(q)q.s.SF_shear,cases));
s0.SF_pin = min(arrayfun(@(q)q.s.SF_pin,cases));
s0.SF_creep = min(arrayfun(@(q)q.s.SF_creep,cases));
s0.SF_wear = min(arrayfun(@(q)q.s.SF_wear,cases));
s0.isValidStress = all(arrayfun(@(q)q.s.isValidStress,cases));

b0.backlash_deg = max(arrayfun(@(q)q.b.backlash_deg,cases));
b0.backlash_rad = b0.backlash_deg*pi/180;
b0.loaded_lost_motion_deg = max(arrayfun(@(q)q.b.loaded_lost_motion_deg,cases));
b0.loaded_lost_motion_rad = b0.loaded_lost_motion_deg*pi/180;
b0.kt_Nm_per_rad = min(arrayfun(@(q)q.b.kt_Nm_per_rad,cases));
b0.isValidBacklash = all(arrayfun(@(q)q.b.isValidBacklash,cases));

reactionNorm = arrayfun(@(q)norm(q.ls.bearingReaction),cases);
[~,iReaction] = max(reactionNorm);
a.s = s0;
a.b = b0;
a.phases = phases;
a.allGeometryValid = all(arrayfun(@(q)q.g.isValidGeometry,cases));
a.allContactsValid = all(arrayfun(@(q)q.ls.isValidContact,cases));
a.minRingContactCount = min(arrayfun(@(q)q.ls.nRingAct,cases));
a.minOutputContactCount = min(arrayfun(@(q)q.ls.nOutAct,cases));
a.maxMinRingGap = max(arrayfun(@(q)q.g.minRingGap,cases));
a.minMinRingGap = min(arrayfun(@(q)q.g.minRingGap,cases));
a.maxTorqueBalanceError = max(arrayfun(@(q) ...
    abs(q.ls.torqueRecovered-p.Tout_nominal),cases));
a.maxBearingReaction_N = cases(iReaction).ls.bearingReaction;
a.cases = cases;
end
