function f = objectiveWrapper(x,Zr,Zc,nOut,p)
%OBJECTIVEWRAPPER Four clean objectives; feasibility stays in constraints.

d = evaluateDesignCached(x,Zr,Zc,nOut,p);
if ~d.isValid || any(~isfinite([d.mass_kg,d.backlash_deg, ...
        d.sigmaContact_MPa,d.kt_Nm_per_rad]))
    f = [1e3,1e3,1e6,1e6];
else
    f = [d.mass_kg,d.backlash_deg,d.sigmaContact_MPa,-d.kt_Nm_per_rad];
end
end
