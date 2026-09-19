function b = backlashModel(g,ls,p)
%BACKLASHMODEL Separate no-load lost motion from elastic torsional stiffness.
% The clearance conversion is a documented analytical surrogate; final
% designs require reverse-contact kinematics and prototype correlation.

if ~g.isValidGeometry || ~ls.isValidContact
    b = invalidBacklash();
    return
end
if strcmpi(p.buildPlane,'XY')
    Edisc = p.E_xy;
else
    Edisc = p.E_z;
end
if p.useMetalPins
    Epin = p.E_pin_metal;
else
    Epin = p.E_pin_printed;
end

wear = p.initialWearAllowance;
ringLever = median(ls.ringMomentArms);
% Geometric clearances are nominal. Process tolerance is sampled separately
% in robustnessModel and is not added a second time here.
ringGap = g.ringClearance+wear;
outputGap = g.outputClearance+wear;
deadbandRing = 2*ringGap/ringLever;
deadbandOutput = 2*outputGap/g.rOutPitch;
deadbandGeom = deadbandRing+deadbandOutput;
thetaCreepAllowance = p.creepBacklashFrac*deadbandGeom;
backlashRad = deadbandGeom+thetaCreepAllowance;

% Elastic compliances in series: ring contact, output-hole material, pins.
kRing = ls.ringTorsionalStiffness;
kLinOut = Edisc*(g.tDisc*g.dOut)/g.minLigament;
kOut = ls.nOutAct*kLinOut*g.rOutPitch^2;

Iring = pi*g.dRing^4/64;
Iout = pi*g.dOut^4/64;
kPinRingLin = 3*Epin*Iring/p.Lpin_ring^3;
kPinOutLin = 3*Epin*Iout/p.Lpin_out^3;
kPin = ls.nRingAct*kPinRingLin*ringLever^2 + ...
       ls.nOutAct*kPinOutLin*g.rOutPitch^2;

stiffnesses = [kRing,kOut,kPin];
if any(~isfinite(stiffnesses)) || any(stiffnesses<=0)
    b = invalidBacklash();
    return
end
kt = 1/sum(1./stiffnesses);
thetaElastic = abs(ls.Tout)/kt;

b.isValidBacklash = isfinite(backlashRad) && isfinite(kt) && kt>0;
b.deadband_ring_rad = deadbandRing;
b.deadband_output_rad = deadbandOutput;
b.deadband_geom_rad = deadbandGeom;
b.creep_allowance_rad = thetaCreepAllowance;
b.backlash_rad = backlashRad;
b.backlash_deg = backlashRad*180/pi;
b.theta_elastic_rad = thetaElastic;
b.loaded_lost_motion_rad = backlashRad+thetaElastic;
b.loaded_lost_motion_deg = b.loaded_lost_motion_rad*180/pi;
b.kt_Nm_per_rad = kt;
b.kRing_Nm_per_rad = kRing;
b.kOutput_Nm_per_rad = kOut;
b.kPins_Nm_per_rad = kPin;
b.modelNote = ['Clearance lost motion is a screening surrogate; ', ...
    'validate using reverse-contact kinematics.'];
end

function b = invalidBacklash()
b.isValidBacklash = false;
b.deadband_ring_rad = inf;
b.deadband_output_rad = inf;
b.deadband_geom_rad = inf;
b.creep_allowance_rad = inf;
b.backlash_rad = inf;
b.backlash_deg = inf;
b.theta_elastic_rad = inf;
b.loaded_lost_motion_rad = inf;
b.loaded_lost_motion_deg = inf;
b.kt_Nm_per_rad = 0;
b.kRing_Nm_per_rad = 0;
b.kOutput_Nm_per_rad = 0;
b.kPins_Nm_per_rad = 0;
b.modelNote = 'Invalid geometry/contact.';
end
