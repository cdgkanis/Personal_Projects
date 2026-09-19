function results = main_optimize()
%MAIN_OPTIMIZE Reproducible mixed-integer multiobjective optimization.

projectDir = fileparts(mfilename('fullpath'));
addpath(projectDir);
if exist('gamultiobj','file')~=2
    error('CycloidalDrive:MissingToolbox', ...
        'gamultiobj is unavailable. Install/enable Global Optimization Toolbox.');
end

p = defaultParams();
rng(p.randomSeed,'twister');
if p.UseParallel
    warning('CycloidalDrive:ParallelReproducibility', ...
        ['Parallel GA population generation is not guaranteed reproducible. ', ...
         'UseParallel=false is recommended for auditable runs.']);
end

% y = [eight continuous variables, index into Zr_set, index into nOut_set].
lb = [p.Rp_min,p.e_min,p.dRing_min,p.dOut_min,p.tDisc_min, ...
      p.clear_min,p.rHole_min,p.outClear_min,1,1];
ub = [p.Rp_max,p.e_max,p.dRing_max,p.dOut_max,p.tDisc_max, ...
      p.clear_max,p.rHole_max,p.outClear_max,numel(p.Zr_set),numel(p.nOut_set)];
intcon = [9,10];
obj = @(y) objectiveMixedWrapper(y,p);
nonlcon = @(y) constraintMixedWrapper(y,p);
[liveOutput,getLiveHistory] = createLiveFeedback(p);
if ~p.liveFeedback
    liveOutput = [];
end

opts = optimoptions('gamultiobj', ...
    'PopulationSize',p.PopulationSize, ...
    'MaxGenerations',p.MaxGenerations, ...
    'CrossoverFcn','crossovertwopoint', ...
    'MutationFcn','mutationadaptfeasible', ...
    'MaxStallGenerations',p.MaxStallGenerations, ...
    'FunctionTolerance',p.FunctionTolerance, ...
    'ConstraintTolerance',p.ConstraintTolerance, ...
    'ParetoFraction',p.ParetoFraction, ...
    'UseParallel',p.UseParallel, ...
    'OutputFcn',liveOutput, ...
    'Display','iter');

[Y,Fsolver,exitflag,output,population,scores] = gamultiobj( ...
    obj,10,[],[],[],[],lb,ub,nonlcon,intcon,opts);

template = struct('x',[],'mixedX',[],'F',[],'design',[]);
allCandidates = repmat(template,size(Y,1),1);
pValidation = p;
pValidation.Ntheta = p.NthetaValidation;
pValidation.nPhaseSamples = p.nPhaseSamplesValidation;
pValidation.cacheToken = p.cacheToken+1;

for i = 1:size(Y,1)
    [x,Zr,Zc,nOut] = decodeMixedDesign(Y(i,:),p);
    d = evaluateDesign(x,Zr,Zc,nOut,pValidation);
    robust = robustnessModel(x,Zr,Zc,nOut,pValidation);
    d.robustFailureProb = robust.failProb;
    d.robustFailureProbLower95 = robust.failProbLower95;
    d.robustFailureProbUpper95 = robust.failProbUpper95;
    d.robustMeanBacklash_deg = robust.meanBacklash_deg;
    d.robustMeanMinSF = robust.meanMinSF;
    d.robustPassByUpperConfidence = robust.passByUpperConfidence;
    d.robustFailureReasons = robust.failureReasons;
    d.robustSampleCount = robust.sampleCount;

    allCandidates(i).x = x;
    allCandidates(i).mixedX = Y(i,:);
    allCandidates(i).F = objectiveWrapper(x,Zr,Zc,nOut,pValidation);
    allCandidates(i).design = d;
end

results = rankParetoSet(allCandidates,pValidation);
results.parameters = pValidation;
results.solver.exitflag = exitflag;
results.solver.output = output;
results.solver.population = population;
results.solver.scores = scores;
results.solver.returnedObjectives = Fsolver;
results.liveHistory = getLiveHistory();
% Exclude the callback closure (which owns graphics handles) from MAT files.
savedOptions = optimoptions(opts,'OutputFcn',[]);
results.manifest = buildRunManifest(projectDir,pValidation,savedOptions,exitflag,output);

resultDir = fullfile(projectDir,'results');
if ~exist(resultDir,'dir')
    mkdir(resultDir);
end
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
resultFile = fullfile(resultDir,['cycloidal_results_',stamp,'.mat']);
results.resultFile = resultFile;
save(resultFile,'results','-v7.3');
fprintf('Saved auditable optimization results to %s\n',resultFile);
end
