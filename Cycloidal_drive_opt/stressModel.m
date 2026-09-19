function s = stressModel(g,ls,p)
%STRESSMODEL Screening stresses with Hertz line contact and pin checks.
% Root bending remains an explicitly reported analytical surrogate and
% should be calibrated against FEA/test data before release.

if ~g.isValidGeometry || ~ls.isValidContact || ...
        g.minLigament<=0 || g.tDisc<=0
    s = invalidStress();
    return
end

t = g.tDisc;
Kload = p.Kover*p.Kdyn*p.Kdist*p.Krel;
Fring = ls.Fring(:)*Kload;
Fout = ls.Fout(:)*Kload;
FrMax = max(Fring);
FoMax = max(Fout);

if strcmpi(p.buildPlane,'XY')
    sigY = p.sigmaY_xy;
    Edisc = p.E_xy;
else
    sigY = p.sigmaY_z;
    Edisc = p.E_z;
end
materialDerate = max(p.Ktemp*p.Korient,realmin);
sigAllowBend = sigY*p.printFactor/materialDerate;
sigAllowBearing = p.sigmaAllow_bearing*p.printFactor/materialDerate;
sigAllowContact = p.sigmaAllow_contact*p.printFactor/materialDerate;
tauAllow = p.tauAllow_shear*p.printFactor/materialDerate;

if p.useMetalPins
    Epin = p.E_pin_metal;
    nuPin = p.nu_pin_metal;
    sigAllowPin = p.sigmaAllow_pin_metal;
else
    Epin = p.E_pin_printed;
    nuPin = p.nu_pin_printed;
    sigAllowPin = p.sigmaAllow_pin_printed*p.printFactor/materialDerate;
end

% Hertz line-contact pressure. The sum of curvatures is conservative when
% local convex/concave classification is uncertain.
rr = g.dRing/2;
rhoProfile = max(ls.ringContactRho(:),p.minCurvatureRadius);
Req = 1./(1/rr+1./rhoProfile);
Econtact = 1/((1-p.nu_petg^2)/Edisc+(1-nuPin^2)/Epin);
sigmaContactEach = sqrt(Fring*Econtact./(pi*t.*Req));
sigmaContact = max(sigmaContactEach);

% Projected bearing pressure at output holes.
sigmaBearingOut = FoMax/(t*g.dOut);

% Root-bending screening surrogate, using actual positive geometry.
rootPitchWidth = g.ringPinSpacing-g.dRing;
wRoot = min(rootPitchWidth,g.lig_hole_foot);
if wRoot<=0
    s = invalidStress();
    return
end
Mroot = FrMax*max(g.e,g.dRing/2);
Iroot = t*wRoot^3/12;
sigmaRoot = Mroot*(wRoot/2)/Iroot;

% Web shear and output-hole ligament equivalent stress.
rEff = g.rFoot;
Aweb = 2*pi*rEff*t;
tauWeb = abs(ls.Tout)/(Aweb*rEff);
wLig = g.minLigament;
sigmaLig = FoMax/(t*wLig);
sigmaVM = sqrt(sigmaLig^2+3*tauWeb^2);
epsRoot = sigmaRoot/Edisc;

% Pin bending is checked in both metal and printed modes.
sigmaPinRing = 32*FrMax*p.Lpin_ring/(pi*g.dRing^3);
sigmaPinOut = 32*FoMax*p.Lpin_out/(pi*g.dOut^3);
sigmaPin = max(sigmaPinRing,sigmaPinOut);

SF_bearing = min(sigAllowContact/sigmaContact, ...
                 sigAllowBearing/sigmaBearingOut);
SF_bend = min(sigAllowBend/sigmaRoot,sigAllowBend/sigmaVM);
SF_shear = tauAllow/tauWeb;
SF_pin = sigAllowPin/sigmaPin;
SF_creep = p.creepStressFrac*sigAllowBend/max(sigmaRoot,sigmaVM);
SF_wear = p.wearStressFrac*min(sigAllowContact/sigmaContact, ...
                               sigAllowBearing/sigmaBearingOut);

s.isValidStress = all(isfinite([sigmaContact,sigmaBearingOut,sigmaRoot, ...
    tauWeb,sigmaLig,sigmaVM,sigmaPin])) && ...
    all([sigmaContact,sigmaBearingOut,sigmaRoot,tauWeb,sigmaLig,sigmaVM,sigmaPin]>=0);
s.sigmaContact_MPa = sigmaContact/1e6;
s.sigmaContactEach_MPa = sigmaContactEach/1e6;
s.sigmaBearingOut_MPa = sigmaBearingOut/1e6;
s.sigmaRoot_MPa = sigmaRoot/1e6;
s.tauWeb_MPa = tauWeb/1e6;
s.sigmaLig_MPa = sigmaLig/1e6;
s.sigmaVM_MPa = sigmaVM/1e6;
s.sigmaPinRing_MPa = sigmaPinRing/1e6;
s.sigmaPinOut_MPa = sigmaPinOut/1e6;
s.epsRoot = epsRoot;
s.effectiveContactRadius_m = min(Req);
s.effectiveContactModulus_Pa = Econtact;
s.rootWidth_m = wRoot;
s.rootBendingModel = 'analytical screening surrogate; calibrate with FEA';
s.SF_bend = SF_bend;
s.SF_bearing = SF_bearing;
s.SF_shear = SF_shear;
s.SF_pin = SF_pin;
s.SF_creep = SF_creep;
s.SF_wear = SF_wear;
end

function s = invalidStress()
fieldsInf = {'sigmaContact_MPa','sigmaBearingOut_MPa','sigmaRoot_MPa', ...
    'tauWeb_MPa','sigmaLig_MPa','sigmaVM_MPa', ...
    'sigmaPinRing_MPa','sigmaPinOut_MPa'};
for i = 1:numel(fieldsInf)
    s.(fieldsInf{i}) = inf;
end
s.sigmaContactEach_MPa = inf;
s.epsRoot = inf;
s.effectiveContactRadius_m = NaN;
s.effectiveContactModulus_Pa = NaN;
s.rootWidth_m = NaN;
s.rootBendingModel = 'invalid geometry/contact';
s.SF_bend = 0;
s.SF_bearing = 0;
s.SF_shear = 0;
s.SF_pin = 0;
s.SF_creep = 0;
s.SF_wear = 0;
s.isValidStress = false;
end
