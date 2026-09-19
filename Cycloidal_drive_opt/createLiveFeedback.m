function [outputFcn,getHistory] = createLiveFeedback(p)
%CREATELIVEFEEDBACK Per-run GA callback with console, figure and saved history.
% No Monte Carlo runs here. Feasibility includes every nominal constraint.
history = struct('generation',{},'elapsedSeconds',{},'feasibleCount',{}, ...
    'populationSize',{},'isFeasible',{},'maxViolation',{},'score',{}, ...
    'mixedX',{},'F',{},'design',{});
started = tic;
fig = [];
geometryAxes = [];
progressAxes = [];
info = [];
graphicsDisabled = ~p.livePlot;
outputFcn = @update;
getHistory = @snapshot;

    function [state,options,optchanged] = update(options,state,flag)
        optchanged = false;
        if ~any(strcmp(flag,{'init','iter','done'})), return; end
        if isempty(state.Population), return; end
        % MATLAB may call 'done' for the last generation a second time.
        if ~isempty(history) && history(end).generation==state.Generation
            return
        end
        n = size(state.Population,1);
        F = inf(n,4);
        violation = inf(n,1);
        designs = cell(n,1);
        % Reuse duplicate members within a generation. This diagnostic pass
        % adds nominal evaluations; it never changes the solver population.
        [members,~,memberMap] = unique(state.Population,'rows','stable');
        for j = 1:size(members,1)
            [x,Zr,Zc,nOut] = decodeMixedDesign(members(j,:),p);
            d = evaluateDesignCached(x,Zr,Zc,nOut,p);
            c = constraintWrapper(x,Zr,Zc,nOut,p);
            rows = find(memberMap==j);
            raw = [d.mass_kg,d.backlash_deg,d.sigmaContact_MPa,-d.kt_Nm_per_rad];
            F(rows,:) = repmat(raw,numel(rows),1);
            violation(rows) = max([0;c(:)]);
            designs(rows) = repmat({d},numel(rows),1);
        end
        choice = selectLiveDesign(F,violation,p.rankWeights,p.ConstraintTolerance);
        if isempty(choice.index), return; end
        i = choice.index;
        d = designs{i};
        entry = struct('generation',state.Generation,'elapsedSeconds',toc(started), ...
            'feasibleCount',choice.feasibleCount,'populationSize',n, ...
            'isFeasible',choice.isFeasible,'maxViolation',choice.maxViolation, ...
            'score',choice.score,'mixedX',state.Population(i,:),'F',F(i,:), ...
            'design',d);
        history(end+1) = entry;
        if choice.isFeasible
            label = 'BEST NOMINAL FEASIBLE';
        else
            label = 'NO FEASIBLE DESIGN: LEAST VIOLATING';
        end
        fprintf(['[Live gen %d] %s | feasible %d/%d | max violation %.3g\n', ...
            '  mass %.4g kg | backlash %.4g deg | contact %.4g MPa | stiffness %.4g Nm/rad\n', ...
            '  Zr/Zc/out %d/%d/%d | Rp/e/dRing/dOut/t/rOut [mm]: %.3f %.3f %.3f %.3f %.3f %.3f\n', ...
            '  ring/output clearance [mm]: %.3f %.3f\n'], ...
            entry.generation,label,entry.feasibleCount,n,entry.maxViolation, ...
            d.mass_kg,d.backlash_deg,d.sigmaContact_MPa,d.kt_Nm_per_rad, ...
            d.Zr,d.Zc,d.nOut,1e3*d.Rp,1e3*d.e,1e3*d.dRing,1e3*d.dOut, ...
            1e3*d.tDisc,1e3*d.rOutPitch,1e3*d.ringClearance,1e3*d.outputClearance);
        if graphicsDisabled, return; end
        % Closing the window disables plotting for the rest of this run.
        if ~isempty(fig) && ~isgraphics(fig)
            graphicsDisabled = true;
            return
        end
        try
            if isempty(fig)
                fig = figure('Name','Cycloidal drive - live generation feedback', ...
                    'NumberTitle','off','Color','w','Position',[80 80 1100 660]);
                geometryAxes = axes('Parent',fig,'Position',[.07 .38 .40 .55]);
                progressAxes = axes('Parent',fig,'Position',[.57 .53 .38 .37]);
                info = uicontrol('Parent',fig,'Style','text','Units','normalized', ...
                    'Position',[.04 .02 .92 .29],'BackgroundColor','w', ...
                    'HorizontalAlignment','left','FontSize',11);
            end
            drawGeometry(geometryAxes,d,p);
            title(geometryAxes,sprintf('Generation %d | %s',entry.generation,label), ...
                'Interpreter','none','FontSize',10);
            generations = [history.generation];
            cla(progressAxes);
            plot(progressAxes,generations,[history.maxViolation],'-o','LineWidth',1.3);
            hold(progressAxes,'on');
            plot(progressAxes,generations, ...
                p.ConstraintTolerance*ones(size(generations)),'--');
            hold(progressAxes,'off');
            xlabel(progressAxes,'Generation');
            ylabel(progressAxes,'Maximum normalized constraint violation');
            title(progressAxes,'Selected candidate in each generation');
            grid(progressAxes,'on');
            set(info,'String',sprintf([ ...
                '%s | generation %d | feasible %d/%d | elapsed %.1f min\n', ...
                'Mass %.4g kg   Backlash %.4g deg   Contact %.4g MPa   Stiffness %.4g Nm/rad\n', ...
                'Zr=%d  Zc=%d  output pins=%d  |  Rp=%.3f  e=%.3f  dRing=%.3f  dOut=%.3f  t=%.3f  rOut=%.3f mm\n', ...
                'Clearances: ring=%.3f mm / output=%.3f mm | min ligament=%.3f mm\n', ...
                'SF: bend %.3g | bearing %.3g | shear %.3g | pin %.3g | creep %.3g | wear %.3g\n', ...
                'Best = weighted compromise on this generation''s nominal feasible Pareto front.\n', ...
                'Robust validation follows optimization. Closing this window keeps console feedback active.'], ...
                label,entry.generation,entry.feasibleCount,n,entry.elapsedSeconds/60, ...
                d.mass_kg,d.backlash_deg,d.sigmaContact_MPa,d.kt_Nm_per_rad, ...
                d.Zr,d.Zc,d.nOut,1e3*d.Rp,1e3*d.e,1e3*d.dRing, ...
                1e3*d.dOut,1e3*d.tDisc,1e3*d.rOutPitch, ...
                1e3*d.ringClearance,1e3*d.outputClearance,1e3*d.minLigament, ...
                d.SF_bend,d.SF_bearing,d.SF_shear,d.SF_pin,d.SF_creep,d.SF_wear));
            drawnow;
        catch ME
            graphicsDisabled = true;
            warning('CycloidalDrive:LivePlotUnavailable', ...
                'Live plot disabled; console/history continue: %s',ME.message);
        end
    end

    function data = snapshot()
        data = history;
    end
end

function drawGeometry(ax,d,p)
g = cycloidGeometry(d.Rp,d.e,d.dRing,d.dOut,d.tDisc, ...
    d.ringClearance,d.rOutPitch,d.Zr,d.Zc,d.nOut,p,d.outputClearance);
cla(ax);
plot(ax,1e3*[g.x;g.x(1)],1e3*[g.y;g.y(1)],'b-','LineWidth',1.4);
hold(ax,'on');
t = linspace(0,2*pi,61);
for j = 1:size(g.ringCtrDisc,1)
    plot(ax,1e3*(g.ringCtrDisc(j,1)+g.dRing/2*cos(t)), ...
        1e3*(g.ringCtrDisc(j,2)+g.dRing/2*sin(t)),'Color',[.25 .25 .25]);
end
for j = 1:size(g.outCtr,1)
    plot(ax,1e3*(g.outCtr(j,1)+g.rHoleGeom*cos(t)), ...
        1e3*(g.outCtr(j,2)+g.rHoleGeom*sin(t)),'Color',[.85 .35 .1]);
end
hold(ax,'off');
axis(ax,'equal');
xlabel(ax,'Disc-frame x (mm)');
ylabel(ax,'Disc-frame y (mm)');
grid(ax,'on');
end
