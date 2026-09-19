function p = defaultParams()
%DEFAULTPARAMS Central parameter set for cycloidal-drive optimization.
% All dimensional values use SI units unless the field name says otherwise.

% Geometry sampling. Increase NthetaValidation for final CAD export/checks.
p.Ntheta = 720;
p.NthetaRobust = 360;
p.NthetaValidation = 2400;
p.inputPhase = 0;
p.nPhaseSamples = 9;
p.nPhaseSamplesRobust = 9;
p.nPhaseSamplesValidation = 31;
p.minDerivativeSpeed = 1e-8;
p.minCurvatureRadius = 0.25e-3;
p.shortWidthCoefficientMin = 0.05;
p.shortWidthCoefficientMax = 0.95;
p.ringPinAssemblyMargin = 0.15e-3;
p.contactEngagementBand = 0.60e-3;
p.maxNoLoadRingGap = 1.00e-3;
p.maxRingInterference = 0.05e-3;
p.minContactMomentArm = 0.05e-3;
p.contactGapScale = 0.15e-3;
p.contactComplianceLength = 2.0e-3;

% PETG material model. Replace with coupon data for the actual process.
p.rho = 1270;
p.E_xy = 2.78e9;
p.E_z  = 2.55e9;
p.nu_petg = 0.38;
p.sigmaY_xy = 75e6;
p.sigmaY_z  = 56e6;
p.sigmaAllow_bearing = 28e6;
p.sigmaAllow_contact = 28e6;
p.tauAllow_shear = 14e6;

% Pin properties are checked for both metal and printed-pin modes.
p.useMetalPins = true;
p.E_pin_metal = 200e9;
p.nu_pin_metal = 0.30;
p.sigmaAllow_pin_metal = 250e6;
p.E_pin_printed = 2.55e9;
p.nu_pin_printed = 0.38;
p.sigmaAllow_pin_printed = 35e6;

% FDM/process factors.
p.printFactor = 0.55;
p.xyTol = 0.15e-3;
p.layerHeight = 0.20e-3;
p.initialWearAllowance = 0.03e-3;
p.buildPlane = 'XY';

% Load multipliers. Temperature/orientation derate allowables only.
p.Kover = 1.25;
p.Kdyn  = 1.15;
p.Kdist = 1.10;
p.Krel  = 1.00;
p.Ktemp = 1.05;
p.Korient = 1.00;

% Long-term screening factors (not service-life predictions).
p.creepStressFrac = 0.35;
p.wearStressFrac = 0.30;
p.creepBacklashFrac = 0.20;

% Minimum safety factors.
p.SFmin_bend = 2.0;
p.SFmin_bearing = 1.5;
p.SFmin_shear = 2.0;
p.SFmin_pin = 2.0;
p.SFmin_creep = 1.3;
p.SFmin_wear = 1.2;

% Geometry/printability.
p.minWall = 2.0e-3;
p.minLigamentReq = 3.0e-3;
p.tDisc_min_req = 4.0e-3;
p.Lpin_ring = 8e-3;
p.Lpin_out = 8e-3;

% Operating point.
p.Tout_nominal = 15; % N*m
p.backlashMax_deg = 1.0;

% Robust validation uses fixed common random numbers and Wilson confidence.
p.robustN = 200;
p.robustSeed = 271828;
p.maxFailureProb = 0.10;
p.robustConfidenceZ = 1.96;
p.torqueSigmaFrac = 0.05;
p.materialSigmaFrac = 0.10;
p.modulusSigmaFrac = 0.08;

% Discrete design sets, optimized as integer indices in one run.
p.Zr_set = 15:2:31;
p.nOut_set = 5:13;

% Continuous bounds: [Rp,e,dRing,dOut,tDisc,ringClr,rOutPitch,outClr].
p.Rp_min = 20e-3;      p.Rp_max = 120e-3;
p.e_min = 0.2e-3;      p.e_max = 2.0e-3;
p.dRing_min = 2e-3;    p.dRing_max = 12e-3;
p.dOut_min = 3e-3;     p.dOut_max = 14e-3;
p.tDisc_min = 4e-3;    p.tDisc_max = 12e-3;
p.clear_min = 0.20e-3; p.clear_max = 0.35e-3;
p.rHole_min = 5e-3;    p.rHole_max = 60e-3;
p.outClear_min = 0.10e-3;
p.outClear_max = 0.35e-3;
p.outputClear_nom = 0.20e-3; % Backward-compatible 7-value x.

% Multiobjective optimizer settings.
p.randomSeed = 314159;
p.PopulationSize = 300;
p.MaxGenerations = 300;
p.MaxStallGenerations = 60;
p.FunctionTolerance = 1e-4;
p.ConstraintTolerance = 1e-6;
p.ParetoFraction = 0.50;
p.UseParallel = true; % Reproducible default.
p.rankWeights = [0.30 0.25 0.25 0.20];
p.liveFeedback = true;
p.livePlot = true; % false keeps console feedback and saved history.
p.cacheToken = 1;

validateParams(p);
end
