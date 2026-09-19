function choice = selectLiveDesign(F,violation,weights,tolerance)
%SELECTLIVEDESIGN Rank feasible nondominated designs in this generation.
% If none is feasible, report the least-violating candidate explicitly.
choice = struct('index',[],'isFeasible',false,'feasibleCount',0, ...
    'score',NaN,'maxViolation',NaN);
if isempty(F), return; end
violation = violation(:);
violation(~isfinite(violation)) = inf;
finite = all(isfinite(F),2);
feasible = finite & violation<=tolerance;
choice.feasibleCount = sum(feasible);
if ~any(feasible)
    [choice.maxViolation,choice.index] = min(violation);
    return
end
idx = find(feasible);
A = F(idx,:);
keep = true(numel(idx),1);
for i = 1:numel(idx)
    keep(i) = ~any(all(A<=A(i,:),2) & any(A<A(i,:),2));
end
idx = idx(keep);
A = F(idx,:);
span = max(A,[],1)-min(A,[],1);
span(span==0) = 1;
normalized = (A-min(A,[],1))./span;
scores = normalized*(weights(:)/sum(weights));
[choice.score,k] = min(scores);
choice.index = idx(k);
choice.isFeasible = true;
choice.maxViolation = violation(choice.index);
end
