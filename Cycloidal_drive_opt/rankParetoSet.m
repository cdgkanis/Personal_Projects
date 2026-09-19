function results = rankParetoSet(allCandidates,p)
%RANKPARETOSET Robust-feasible global nondominated filtering and ranking.

results.allCandidates = allCandidates;
results.allPareto = struct('x',{},'F',{},'design',{});
results.score = [];
results.sortedIdx = [];
results.best = [];
results.status = 'no candidates';
if isempty(allCandidates)
    return
end

n = numel(allCandidates);
F = NaN(n,4);
feasible = false(n,1);
for i = 1:n
    F(i,:) = allCandidates(i).F;
    d = allCandidates(i).design;
    feasible(i) = d.isValid && all(isfinite(F(i,:))) && ...
        isfinite(d.robustFailureProbUpper95) && ...
        d.robustFailureProbUpper95<=p.maxFailureProb;
end
candidateIdx = find(feasible);
if isempty(candidateIdx)
    results.status = 'no robust-feasible candidates';
    results.feasibleMask = feasible;
    return
end

Ffeasible = F(candidateIdx,:);
keep = true(size(candidateIdx));
for i = 1:numel(candidateIdx)
    dominatesI = all(Ffeasible<=Ffeasible(i,:),2) & ...
                 any(Ffeasible<Ffeasible(i,:),2);
    dominatesI(i) = false;
    if any(dominatesI)
        keep(i) = false;
    end
end
paretoIdx = candidateIdx(keep);
front = allCandidates(paretoIdx);
Fp = F(paretoIdx,:);

fmin = min(Fp,[],1);
fmax = max(Fp,[],1);
span = fmax-fmin;
span(span<=sqrt(eps)) = 1;
Fn = (Fp-fmin)./span;
w = p.rankWeights(:)/sum(p.rankWeights);
score = Fn*w;
[~,order] = sort(score,'ascend');

results.allPareto = front;
results.globalParetoOriginalIdx = paretoIdx;
results.feasibleMask = feasible;
results.normalizedObjectives = Fn;
results.score = score;
results.sortedIdx = order;
results.best = front(order(1));
results.status = 'ok';
results.rankWeights = w.';
end
