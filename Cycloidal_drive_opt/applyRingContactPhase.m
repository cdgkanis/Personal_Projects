function gp = applyRingContactPhase(g,p,alpha,Zc)
%APPLYRINGCONTACTPHASE Reassemble ring pins in the disc frame at one phase.

discAngle = -alpha/Zc;
R = [cos(discAngle),-sin(discAngle);sin(discAngle),cos(discAngle)];
discCentreHousing = g.e*[cos(alpha),sin(alpha)];
ringCtrDisc = (g.ringCtr-discCentreHousing)*R;
Zr = size(ringCtrDisc,1);
rr = g.dRing/2;

pinGap = zeros(Zr,1);
idxClosest = zeros(Zr,1);
forceNormal = zeros(Zr,2);
momentArm = zeros(Zr,1);
profileAngle = atan2(g.y,g.x);
for i = 1:Zr
    pinAngle = atan2(ringCtrDisc(i,2),ringCtrDisc(i,1));
    angleDelta = abs(atan2(sin(profileAngle-pinAngle), ...
                           cos(profileAngle-pinAngle)));
    searchIdx = find(angleDelta<=2*pi/Zc);
    if isempty(searchIdx), searchIdx = (1:numel(g.x)).'; end
    delta = [g.x(searchIdx)-ringCtrDisc(i,1), ...
             g.y(searchIdx)-ringCtrDisc(i,2)];
    dist = hypot(delta(:,1),delta(:,2));
    [centreDistance,localIdx] = min(dist);
    idx = searchIdx(localIdx);
    pinGap(i) = centreDistance-rr;
    idxClosest(i) = idx;
    if centreDistance>p.minDerivativeSpeed
        normal = [g.x(idx)-ringCtrDisc(i,1), ...
                  g.y(idx)-ringCtrDisc(i,2)]/centreDistance;
    else
        normal = [NaN,NaN];
    end
    forceNormal(i,:) = normal;
    momentArm(i) = g.x(idx)*normal(2)-g.y(idx)*normal(1);
end

gp = g;
gp.inputPhase = alpha;
gp.ringCtrDisc = ringCtrDisc;
gp.pinGap = pinGap;
gp.pinClosestProfileIdx = idxClosest;
gp.pinForceNormal = forceNormal;
gp.pinMomentArm = momentArm;
gp.minRingGap = min(pinGap);
gp.contactReachable = gp.minRingGap<=p.maxNoLoadRingGap && ...
    gp.minRingGap>=-p.maxRingInterference;
gp.isValidGeometry = g.isValidGeometryBase && gp.contactReachable;
end
