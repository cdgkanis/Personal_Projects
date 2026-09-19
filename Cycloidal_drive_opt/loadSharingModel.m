function ls = loadSharingModel(g,p)
%LOADSHARINGMODEL Deterministic phase-dependent contact/load distribution.
% Ring loads are nonnegative normal forces and exactly balance output torque.
% The eccentric bearing carries the residual vector force.

T = p.Tout_nominal;
torqueSign = sign(T);
if torqueSign==0
    torqueSign = 1;
end

gap = g.pinGap(:);
leverSigned = g.pinMomentArm(:);
gapExcess = gap-min(gap);
candidate = isfinite(gap) & isfinite(leverSigned) & ...
    gapExcess<=p.contactEngagementBand & ...
    torqueSign*leverSigned>p.minContactMomentArm;
active = find(candidate);

validRingContact = ~isempty(active) && g.contactReachable;
if validRingContact
    lever = abs(leverSigned(active));
    raw = exp(-gapExcess(active)/max(p.contactGapScale,realmin)).*lever;
    if ~all(isfinite(raw)) || sum(raw)<=0
        validRingContact = false;
    end
end

if validRingContact
    share = raw/sum(raw);
    scale = abs(T)/sum(share.*lever);
    Fring = scale*share;
    normals = g.pinForceNormal(active,:);
    ringForceVectors = normals.*Fring;
    torqueRecovered = sum(Fring.*lever)*torqueSign;

    if strcmpi(p.buildPlane,'XY')
        Edisc = p.E_xy;
    else
        Edisc = p.E_z;
    end
    kNormal = Edisc*g.tDisc*g.dRing/max(p.contactComplianceLength,realmin);
    ringDeflection = Fring/max(kNormal,realmin);
    ringTorsionalStiffness = sum(kNormal*lever.^2);
    contactRho = g.rho(g.pinClosestProfileIdx(active));
else
    active = zeros(0,1);
    Fring = zeros(0,1);
    ringForceVectors = zeros(0,2);
    torqueRecovered = 0;
    ringDeflection = zeros(0,1);
    ringTorsionalStiffness = 0;
    contactRho = zeros(0,1);
end

% Output-pin forces are tangential and exactly balance the same torque.
nOut = size(g.outCtr,1);
validOutputContact = nOut>=1 && g.rOutPitch>0 && g.minLigament>0;
if validOutputContact
    % Equal-compliance carrier pins all participate; the small deterministic
    % harmonic represents a documented assembly/load-distribution bias.
    outputActive = (1:nOut).';
    phase = (outputActive-1)*2*pi/nOut;
    phaseOffset = p.inputPhase;
    if isfield(g,'inputPhase')
        phaseOffset = g.inputPhase;
    end
    maldistribution = 1+0.10*cos(phase-phaseOffset);
    outputShare = maldistribution/sum(maldistribution);
    FoTotal = abs(T)/g.rOutPitch;
    Fout = FoTotal*outputShare;
    tangential = torqueSign*[-sin(phase),cos(phase)];
    outputForceVectors = tangential.*Fout;
else
    outputActive = zeros(0,1);
    Fout = zeros(0,1);
    outputForceVectors = zeros(0,2);
end

ls.Tout = T;
ls.isValidContact = validRingContact && validOutputContact;
ls.activeRingPins = active;
ls.activeOutputPins = outputActive;
ls.nRingAct = numel(active);
ls.nOutAct = numel(outputActive);
ls.gap_j = gap(active);
ls.gapExcess = gapExcess(active);
ls.ringMomentArms = abs(leverSigned(active));
ls.ringContactRho = contactRho;
ls.Fring = Fring(:);
ls.Fout = Fout(:);
ls.ringForceVectors = ringForceVectors;
ls.outputForceVectors = outputForceVectors;
ls.bearingReaction = -sum(ringForceVectors,1);
ls.ringDeflection = ringDeflection;
ls.ringTorsionalStiffness = ringTorsionalStiffness;
ls.torqueRecovered = torqueRecovered;

if isempty(Fring), ls.FrMax = inf; else, ls.FrMax = max(Fring); end
if isempty(Fout),  ls.FoMax = inf; else, ls.FoMax = max(Fout); end
end
