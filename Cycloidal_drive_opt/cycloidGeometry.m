function g = cycloidGeometry(Rp,e,dRing,dOut,tDisc,ringClr,rOutPitch,Zr,Zc,nOut,p,varargin)
%CYCLOIDGEOMETRY Conjugate roller-pin cycloidal profile and validity margins.
% Coordinates are in the cycloidal-disc frame. The profile follows the
% standard conjugate pin-gear form with Zr pins and Zc=Zr-1 lobes.

validateattributes([Rp,e,dRing,dOut,tDisc,ringClr,rOutPitch],{'numeric'}, ...
    {'real','finite','positive'});
validateattributes([Zr,Zc,nOut],{'numeric'}, ...
    {'real','finite','integer','positive'});
if Zc~=Zr-1
    error('CycloidalDrive:ToothCountMismatch','Zc must equal Zr-1.');
end
if isempty(varargin)
    outClr = p.outputClear_nom;
else
    outClr = varargin{1};
    validateattributes(outClr,{'numeric'},{'real','finite','positive'});
end

Npts = round(p.Ntheta);
if Npts < 8*Zc
    error('CycloidalDrive:InsufficientSampling', ...
        'Ntheta must provide at least eight samples per cycloidal lobe.');
end
phi = (0:Npts-1)'*(2*pi/Npts);
h = 2*pi/Npts;

rr = dRing/2;
k1 = e*Zr/Rp;
shortWidthMargin = min(k1-p.shortWidthCoefficientMin, ...
    p.shortWidthCoefficientMax-k1);
S = sqrt(max(1 + k1^2 - 2*k1*cos(Zc*phi),realmin));
q = rr + ringClr;

% Conjugate roller-pin profile (a constant normal modification q).
x = (Rp - q./S).*sin(phi) - (e - k1*q./S).*sin(Zr*phi);
y = (Rp - q./S).*cos(phi) - (e - k1*q./S).*cos(Zr*phi);

% Periodic finite differences of the final working profile.
dx  = (circshift(x,-1)-circshift(x,1))/(2*h);
dy  = (circshift(y,-1)-circshift(y,1))/(2*h);
ddx = (circshift(x,-1)-2*x+circshift(x,1))/(h^2);
ddy = (circshift(y,-1)-2*y+circshift(y,1))/(h^2);
spd = hypot(dx,dy);
kappa = (dx.*ddy-dy.*ddx)./max(spd.^3,realmin);
rho = 1./max(abs(kappa),realmin);
finiteRho = rho(isfinite(rho));
if isempty(finiteRho)
    minRho = 0;
else
    minRho = min(finiteRho);
end
minSpeed = min(spd);

rProf = hypot(x,y);
rFoot = min(rProf);

% Verify the expected Zc-fold radial symmetry numerically.
phi2 = [phi; 2*pi];
r2 = [rProf; rProf(1)];
rShift = interp1(phi2,r2,mod(phi+2*pi/Zc,2*pi),'linear');
symmetryError = max(abs(rShift-rProf));

% Use MATLAB's robust polygon representation without silently repairing it.
profileIntersect = 1;
Apoly = NaN;
try
    warningState = warning;
    warningCleanup = onCleanup(@() warning(warningState));
    warning('off','all');
    pg = polyshape(x,y,'Simplify',false,'KeepCollinearPoints',true);
    if issimplified(pg)
        profileIntersect = -1;
        Apoly = area(pg);
    end
    clear warningCleanup
catch
    if exist('warningCleanup','var'), clear warningCleanup; end
    profileIntersect = 1;
end

% Ring pin layout and overlap margin.
angRing = (0:Zr-1)'*2*pi/Zr;
ringCtrHousing = [Rp*cos(angRing), Rp*sin(angRing)];
ringPinSpacing = 2*Rp*sin(pi/Zr);
ringPinSpacingMargin = ringPinSpacing-dRing-p.ringPinAssemblyMargin;

% Output holes. rOutPitch is their centre-circle radius.
angOut = (0:nOut-1)'*2*pi/nOut;
outCtr = [rOutPitch*cos(angOut),rOutPitch*sin(angOut)];
dHole = dOut + 2*e + 2*outClr;
rHoleGeom = dHole/2;
rHoleOuter = rOutPitch+rHoleGeom;
outputHoleViolation = rHoleOuter-rFoot;
if nOut>1
    holeCenterDist = 2*rOutPitch*sin(pi/nOut);
else
    holeCenterDist = inf;
end
lig_hole_hole = holeCenterDist-2*rHoleGeom;
lig_hole_foot = rFoot-rHoleOuter;
minLigament = min(lig_hole_hole,lig_hole_foot);

% Assemble the ring pins in the disc frame at the requested input phase.
alpha = p.inputPhase;
discAngle = -alpha/Zc;
R = [cos(discAngle),-sin(discAngle);sin(discAngle),cos(discAngle)];
discCentreHousing = e*[cos(alpha),sin(alpha)];
ringCtrDisc = (ringCtrHousing-discCentreHousing)*R;

pinGap = zeros(Zr,1);
pinClosestProfileIdx = zeros(Zr,1);
pinForceNormal = zeros(Zr,2);
pinMomentArm = zeros(Zr,1);
profileAngle = atan2(y,x);
for i = 1:Zr
    pinAngle = atan2(ringCtrDisc(i,2),ringCtrDisc(i,1));
    angleDelta = abs(atan2(sin(profileAngle-pinAngle), ...
                           cos(profileAngle-pinAngle)));
    searchIdx = find(angleDelta<=2*pi/Zc);
    if isempty(searchIdx), searchIdx = (1:Npts).'; end
    delta = [x(searchIdx)-ringCtrDisc(i,1), ...
             y(searchIdx)-ringCtrDisc(i,2)];
    dist = hypot(delta(:,1),delta(:,2));
    [centreDistance,localIdx] = min(dist);
    idx = searchIdx(localIdx);
    pinGap(i) = centreDistance-rr;
    pinClosestProfileIdx(i) = idx;
    if centreDistance>p.minDerivativeSpeed
        nForce = [x(idx)-ringCtrDisc(i,1), ...
                  y(idx)-ringCtrDisc(i,2)]/centreDistance;
    else
        nForce = [NaN,NaN];
    end
    pinForceNormal(i,:) = nForce;
    pinMomentArm(i) = x(idx)*nForce(2)-y(idx)*nForce(1);
end
minRingGap = min(pinGap);

geometryFinite = all(isfinite([x;y;kappa])) && isfinite(Apoly);
simpleProfile = profileIntersect<0;
derivativeValid = minSpeed>=p.minDerivativeSpeed;
curvatureMargin = minRho-p.minCurvatureRadius;
symmetryMargin = 1e-3*Rp-symmetryError;
holesContained = outputHoleViolation<=0 && minLigament>=p.minLigamentReq;
pinsSeparated = ringPinSpacingMargin>=0;
contactReachable = minRingGap<=p.maxNoLoadRingGap && ...
    minRingGap>=-p.maxRingInterference;
shortWidthValid = shortWidthMargin>=0;
isValidGeometryBase = geometryFinite && simpleProfile && derivativeValid && ...
    curvatureMargin>=0 && symmetryMargin>=0 && holesContained && ...
    pinsSeparated && shortWidthValid;
isValidGeometry = isValidGeometryBase && contactReachable;

if geometryFinite && simpleProfile && outputHoleViolation<=0 && ...
        lig_hole_hole>=0 && Apoly>0
    Aholes = nOut*pi*rHoleGeom^2;
    Adisc = Apoly-Aholes;
    if Adisc<=0
        Adisc = NaN;
        mass_kg = inf;
    else
        mass_kg = Adisc*tDisc*p.rho;
    end
else
    Aholes = NaN;
    Adisc = NaN;
    mass_kg = inf;
end

g.th = phi;
g.x = x; g.y = y;
g.dx = dx; g.dy = dy;
g.ddx = ddx; g.ddy = ddy;
g.kappa = kappa; g.rho = rho;
g.minCurvatureRadius = minRho;
g.curvatureMargin = curvatureMargin;
g.minDerivativeSpeed = minSpeed;
g.profileSymmetryError = symmetryError;
g.symmetryMargin = symmetryMargin;
g.rProf = rProf;
g.rMin = min(rProf);
g.rMax = max(rProf);
g.rFoot = rFoot;
g.shortWidthCoefficient = k1;
g.shortWidthMargin = shortWidthMargin;

g.Rp = Rp; g.e = e;
g.dRing = dRing; g.dOut = dOut;
g.tDisc = tDisc;
g.clearance = ringClr;
g.ringClearance = ringClr;
g.outputClearance = outClr;
g.rHole = rOutPitch; % Legacy field name.
g.rOutPitch = rOutPitch;
g.dHole = dHole;
g.rHoleGeom = rHoleGeom;
g.ringCtr = ringCtrHousing;
g.ringCtrDisc = ringCtrDisc;
g.outCtr = outCtr;
g.pinGap = pinGap;
g.pinClosestProfileIdx = pinClosestProfileIdx;
g.pinForceNormal = pinForceNormal;
g.pinMomentArm = pinMomentArm;
g.minRingGap = minRingGap;

g.ringPinSpacing = ringPinSpacing;
g.ringPinSpacingMargin = ringPinSpacingMargin;
g.minLigament = minLigament;
g.lig_hole_hole = lig_hole_hole;
g.lig_hole_foot = lig_hole_foot;
g.outputHoleViolation = outputHoleViolation;
g.footViolation = outputHoleViolation; % Legacy alias only.
g.profileIntersect = profileIntersect;
g.contactReachable = contactReachable;
g.isValidGeometryBase = isValidGeometryBase;
g.isValidGeometry = isValidGeometry;
g.areaOuter_m2 = Apoly;
g.areaHoles_m2 = Aholes;
g.area_m2 = Adisc;
g.mass_kg = mass_kg;
end
