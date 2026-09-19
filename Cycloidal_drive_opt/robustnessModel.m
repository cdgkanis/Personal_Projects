function r = robustnessModel(x,Zr,Zc,nOut,p)
%ROBUSTNESSMODEL Fixed-sample post-optimization uncertainty validation.
% Common random numbers make repeated evaluations exactly reproducible.

v = unpackDesign(x,p);
x0 = v.x;
N = round(p.robustN);
stream = RandStream('mt19937ar','Seed',p.robustSeed);
Z = randn(stream,N,12);

fails = 0;
backlashValues = NaN(N,1);
minSFValues = NaN(N,1);
reason.geometry = 0;
reason.contact = 0;
reason.strength = 0;
reason.backlash = 0;
reason.invalidSample = 0;

for i = 1:N
    xp = x0;
    xp(1) = xp(1)+p.xyTol*Z(i,1);
    xp(2) = xp(2)+0.10*p.xyTol*Z(i,2);
    xp(3) = xp(3)+p.xyTol*Z(i,3);
    xp(4) = xp(4)+p.xyTol*Z(i,4);
    xp(5) = xp(5)+p.layerHeight*Z(i,5);
    xp(6) = xp(6)+0.50*p.xyTol*Z(i,6);
    xp(7) = xp(7)+p.xyTol*Z(i,7);
    xp(8) = xp(8)+0.50*p.xyTol*Z(i,8);

    if any(~isfinite(xp)) || any(xp<=0)
        fails = fails+1;
        reason.invalidSample = reason.invalidSample+1;
        continue
    end

    p2 = p;
    p2.Ntheta = p.NthetaRobust;
    p2.nPhaseSamples = p.nPhaseSamplesRobust;
    strengthFactor = exp(p.materialSigmaFrac*Z(i,9)-0.5*p.materialSigmaFrac^2);
    modulusFactor = exp(p.modulusSigmaFrac*Z(i,10)-0.5*p.modulusSigmaFrac^2);
    bearingFactor = exp(0.08*Z(i,11)-0.5*0.08^2);
    torqueFactor = 1+p.torqueSigmaFrac*Z(i,12);
    if torqueFactor<=0
        fails = fails+1;
        reason.invalidSample = reason.invalidSample+1;
        continue
    end
    p2.sigmaY_xy = p.sigmaY_xy*strengthFactor;
    p2.sigmaY_z = p.sigmaY_z*strengthFactor;
    p2.sigmaAllow_bearing = p.sigmaAllow_bearing*bearingFactor;
    p2.sigmaAllow_contact = p.sigmaAllow_contact*bearingFactor;
    p2.tauAllow_shear = p.tauAllow_shear*bearingFactor;
    p2.E_xy = p.E_xy*modulusFactor;
    p2.E_z = p.E_z*modulusFactor;
    p2.Tout_nominal = p.Tout_nominal*torqueFactor;

    try
        ds = evaluateDesign(xp,Zr,Zc,nOut,p2);
    catch
        fails = fails+1;
        reason.invalidSample = reason.invalidSample+1;
        continue
    end

    if ~ds.isValidGeometry
        fails = fails+1;
        reason.geometry = reason.geometry+1;
        continue
    elseif ~ds.isValidContact
        fails = fails+1;
        reason.contact = reason.contact+1;
        continue
    end

    minSF = min([ds.SF_bend,ds.SF_bearing,ds.SF_shear,ds.SF_pin, ...
                 ds.SF_creep,ds.SF_wear]);
    backlashValues(i) = ds.backlash_deg;
    minSFValues(i) = minSF;
    badStrength = ~ds.isValidStress || ...
        ds.SF_bend<p2.SFmin_bend || ...
        ds.SF_bearing<p2.SFmin_bearing || ...
        ds.SF_shear<p2.SFmin_shear || ...
        ds.SF_pin<p2.SFmin_pin || ...
        ds.SF_creep<p2.SFmin_creep || ...
        ds.SF_wear<p2.SFmin_wear;
    badBacklash = ~ds.isValidBacklash || ds.backlash_deg>p2.backlashMax_deg;
    if badStrength || badBacklash
        fails = fails+1;
        reason.strength = reason.strength+double(badStrength);
        reason.backlash = reason.backlash+double(badBacklash);
    end
end

phat = fails/N;
z = p.robustConfidenceZ;
den = 1+z^2/N;
centre = (phat+z^2/(2*N))/den;
half = z*sqrt(phat*(1-phat)/N+z^2/(4*N^2))/den;

r.failCount = fails;
r.sampleCount = N;
r.failProb = phat;
r.failProbLower95 = max(0,centre-half);
r.failProbUpper95 = min(1,centre+half);
r.passByUpperConfidence = r.failProbUpper95<=p.maxFailureProb;
r.meanBacklash_deg = mean(backlashValues,'omitnan');
r.meanMinSF = mean(minSFValues,'omitnan');
r.validBacklashSamples = sum(isfinite(backlashValues));
r.validSafetyFactorSamples = sum(isfinite(minSFValues));
r.failureReasons = reason;
r.seed = p.robustSeed;
r.commonRandomNumbers = true;
end
