classdef regressorManager < matlab.mixin.Copyable
%classdef regressorManager < handle
    properties
        % Genetic algorithm properties
        maxPopulation  % Maximum allowed population size.
        crossoverrate  % Crossover rate.
        mutationrate   % Mutation rate for genetic diversity.
        % Generation estructure:
        regPopulation  % Regressors vector
        x              % Input data.
        y              % Target output values.
        n              % Time index.
        Qpmax          % Maximum noncausal memory.
        Qnmax          % Maximum causal memory.
        Pmax           % Maximum nonlinear order.
        X              % Regressors matrix.
        yX             % Output in sync with U.
        s              % Support set.
        Rmat           % Regressors matrix in text format.
        nmse           % Normalized Mean Squared Error metric of BIC.
        nopt           % Optimum number of regressors (BIC).
        nmsev          % NMSEv
        % Output config:
        verbosity      % Level of output detail.
        showPlots      % Toggle for plotting results.
        inittype       % itilization type
        evaluationtype
    end
    %%
    methods
        %% regressorManager: builder
        function self = regressorManager(x,y,GVGconfig)
            self.maxPopulation = GVGconfig.maxPopulation;
            self.x = x;
            self.y = y;
            self.Qpmax = GVGconfig.Qpmax;
            self.Qnmax = GVGconfig.Qnmax;
            self.Pmax = GVGconfig.Pmax;
            self.n = (1:length(x))';
            self.regPopulation = [];
            self.verbosity = GVGconfig.verbosity;
            self.showPlots = GVGconfig.showPlots;
            self.crossoverrate = GVGconfig.crossoverrate;
            self.mutationrate = GVGconfig.mutationrate;
            self.inittype = GVGconfig.inittype;
            self.evaluationtype = GVGconfig.evaluationtype;
        end
        %% initialization: creates initial population
        %       Create the population with a defined initial set:
        function self = initialization(self)
            modelsConfig;

            % Default: 3 regressors: R([0],[],[]),R([],[0],[]),R([],[],[0])
            if strcmp(self.inittype,'default') || strcmp(self.inittype, 'compositeall')
                self.regPopulation = [self.regPopulation Regressor([0],[],[])];
                self.regPopulation = [self.regPopulation Regressor([],[0],[])];
                self.regPopulation = [self.regPopulation Regressor([],[],[0])];
            % Full Volterra (Volterra series up to 5th order, polinomic up
            % to 13th order).
            end
            
            if strcmp(self.inittype,'FV') || strcmp(self.inittype, 'compositeall')
                reg = fv(Pfv,Mfv);
                for ireg = 1:length(reg.q)
                    X = [];
                    Xconj = [];
                    for iconst = 1:length(reg.q{ireg})
                        if(reg.c{ireg}(iconst))
                            Xconj = [Xconj reg.q{ireg}(iconst)];
                        else
                            X = [X reg.q{ireg}(iconst)];
                        end
                    end
                    self.regPopulation = [self.regPopulation Regressor(X,Xconj,[])];
                    self.regPopulation(ireg).deriveEnvelopeTerms();
                    self.regPopulation(ireg).sortindexes();
                end
            end
            % Complex valued (13th order).
            if strcmp(self.inittype,'CVS') || strcmp(self.inittype, 'compositeall')
                reg = cvs(Pcvs,Mcvs);
                for ireg = 1:length(reg.q)
                    X = [];
                    Xconj = [];
                    for iconst = 1:length(reg.q{ireg})
                        if(reg.c{ireg}(iconst))
                            Xconj = [Xconj reg.q{ireg}(iconst)];
                        else
                            X = [X reg.q{ireg}(iconst)];
                        end
                    end
                    self.regPopulation = [self.regPopulation Regressor(X,Xconj,[])];
                    self.regPopulation(ireg).deriveEnvelopeTerms();
                    self.regPopulation(ireg).sortindexes();
                end
            end
            % Memory polynomial (13th order, memory depth 10).
            if strcmp(self.inittype,'MP') || strcmp(self.inittype, 'compositeall')
                P = Pmp;                
                M = Mmp;
                for k = 0:((P-1)/2)
                    for l = 0:M
                        X = l;
                        Xenv = repmat(l,1,2*k);

                        self.regPopulation = [self.regPopulation Regressor(X,[],Xenv)];
                        self.regPopulation(end).deriveEnvelopeTerms();
                        self.regPopulation(end).sortindexes();
                    end
                end
            end
            % DDR (2nd order dynamics, 13th order, memory depth 10).
            if strcmp(self.inittype,'DDR') || strcmp(self.inittype, 'compositeall')
                % 1st
                P = Pddr;                
                M = Mddr;

                for k = 0:((P-1)/2)
                    for l = 0:M
                        X = l;
                        Xenv = repmat(0,1,2*k);

                        self.regPopulation = [self.regPopulation Regressor(X,[],Xenv)];
                        self.regPopulation(end).deriveEnvelopeTerms();
                        self.regPopulation(end).sortindexes();
                    end
                end
                % 2nd
                for k = 1:((P-1)/2)
                    for l = 1:M
                        X = [0 0];
                        Xconj = l;
                        Xenv = repmat(0,1,2*(k-1));
                        self.regPopulation = [self.regPopulation Regressor(X,Xconj,Xenv)];
                        self.regPopulation(end).deriveEnvelopeTerms();
                        self.regPopulation(end).sortindexes();
                    end
                end
                % 3rd
                for k = 1:((P-1)/2)
                    for l1 = 1:M
                        for l2 = l1:M
                            X = [l1 l2];
                            Xconj = 0;
                            Xenv = repmat(0,1,2*(k-1));
                            self.regPopulation = [self.regPopulation Regressor(X,Xconj,Xenv)];
                            self.regPopulation(end).deriveEnvelopeTerms();
                            self.regPopulation(end).sortindexes();
                        end
                    end
                end
                for k = 1:((P-1)/2)
                    for l1 = 1:M
                        for l2 = 1:M
                            X = [0 l2];
                            Xconj = [l1];
                            Xenv = repmat(0,1,2*(k-1));
                            self.regPopulation = [self.regPopulation Regressor(X,Xconj,Xenv)];
                            self.regPopulation(end).deriveEnvelopeTerms();
                            self.regPopulation(end).sortindexes();
                        end
                    end
                end
                for k = 2:((P-1)/2)
                    for l1 = 1:M
                        for l2 = l1:M
                            X = [0 0 0];
                            Xconj = [l1 l2];
                            Xenv = repmat(0,1,2*(k-2));
                            self.regPopulation = [self.regPopulation Regressor(X,Xconj,Xenv)];
                            self.regPopulation(end).deriveEnvelopeTerms();
                            self.regPopulation(end).sortindexes();
                        end
                    end
                end
            end
            % GMP
            if strcmp(self.inittype,'GMP') || strcmp(self.inittype, 'compositeall')
                PgmpEff = max(1, min(Pgmp, self.Pmax));
                KaEff = 0:(PgmpEff-1);
                KbcEff = 1:(PgmpEff-1);

                for k = 1:length(KaEff)
                    for l = 0:min(Lgmp, self.Qpmax)
                        X = l;
                        Xenv = repmat(l,1,KaEff(k));

                        self.regPopulation = [self.regPopulation Regressor(X,[],Xenv)];
                        self.regPopulation(end).deriveEnvelopeTerms();
                        self.regPopulation(end).sortindexes();
                    end
                end

                for k = 1:length(KbcEff)
                    for l = 0:min(Lgmp, self.Qpmax)
                        mMax = min(Mgmp, self.Qpmax - l);
                        for m = 1:mMax

                            X = l;
                            Xenv = repmat(l+m,1,KbcEff(k));

                            self.regPopulation = [self.regPopulation Regressor(X,[],Xenv)];
                            self.regPopulation(end).deriveEnvelopeTerms();
                            self.regPopulation(end).sortindexes();
                        end
                    end
                end

                for k = 1:length(KbcEff)
                    for l = 0:min(Lgmp, self.Qpmax)
                        mMax = min(Mgmp, l + self.Qnmax);
                        for m = 1:mMax
                            X = l;
                            Xenv = repmat(l-m,1,KbcEff(k));

                            self.regPopulation = [self.regPopulation Regressor(X,[],Xenv)];
                            self.regPopulation(end).deriveEnvelopeTerms();
                            self.regPopulation(end).sortindexes();
                        end
                    end
                end
            end
        end
        %% buildU
        %   Builds the regressor matrix
        function self = buildX(self)
            self.X = [];
            nRows = numel(self.n(1+self.Qpmax:end-self.Qnmax));
            nRegs = length(self.regPopulation);
            estimatedGB = nRows * nRegs * 16 / 2^30;
            if estimatedGB > 0.5
                fprintf('[regressorManager.buildX] Materializando X: rows=%d regs=%d estimado=%.2f GB\n', ...
                    nRows, nRegs, estimatedGB);
            end
            for i = 1:length(self.regPopulation)
                self.regPopulation(i).buildRegressor(self.x, self.n, self.Qpmax, self.Qnmax);
                self.X =  [self.X, self.regPopulation(i).reg];
                self.Rmat{i} = self.regPopulation(i).print();
            end
            self.yX = self.y(self.n(1+self.Qpmax:end-self.Qnmax));
        end
        %% printModel
        %   Generates Rmat in human-readable form
        function self = printModel(self)
            for i = 1:length(self.regPopulation)
                self.Rmat{i} = self.regPopulation(i).print();
                fprintf('%d %s\n',i,self.Rmat{i});
            end
        end

        %% evaluation
        %   Executes domp over the population
        function self = evaluation(self)
            self.buildX();

            [~, s, nopt, ~, ~, nmse] = RCDOMP_GVG(self.X, self.yX, self.Rmat, self.maxPopulation, self.verbosity, self.showPlots, self.evaluationtype);
            self.nmsev = nmse;
            self.nopt = nopt;

            if strcmp(self.evaluationtype,'BIC')
            self.nmse = nmse(nopt);
            self.s = s(1:nopt);
            self.regPopulation = self.regPopulation(s(1:nopt));
            elseif strcmp(self.evaluationtype,'maxPopulation')
                self.nmse = nmse(end);
                self.s = s;
                self.regPopulation = self.regPopulation(s);
            end
            scores = [0 diff(nmse)];

            for i = 1:length(self.regPopulation)
                self.regPopulation(i).score = scores(i);
            end
        end
        %% selection:
        function self = selection(self)
            self.regPopulation = self.regPopulation(1:min(self.maxPopulation,length(self.regPopulation)));
        end
        %% crossover: creates the next generation
        %       Mixes the best regressors between them to create
        %       the next generation of regressors.
        function self = crossover(self)
            rp = randperm(length(self.regPopulation));
            for i = 1:(length(self.regPopulation)*self.crossoverrate)
                newr = self.regPopulation(i).crossover(self.regPopulation(rp(i)));
                if(self.verbosity>=2)
                    fprintf("Crossover: %s and %s produced %s\n", self.regPopulation(i).print(), self.regPopulation(rp(i)).print(), newr.print());
                end
                self.regPopulation = [self.regPopulation newr];
            end
        end
        %% mutation: Mutates part of the populaton,
        %       Mutates part of the populaton,
        %       determinated by the input mutationrate.
        function self = mutation(self)
            Xmutate = self.regPopulation(1:floor(length(self.regPopulation)*self.mutationrate));
            for i = 1:length(Xmutate)
                r = randi([1 3],1,1);
                muttype = ["functional", "memory", "order"];
                if(self.verbosity>=2) fprintf("Mutation (type %s): %s mutated to", muttype(r), Xmutate(i).print()); end
                if r == 1
                    % Functional mutation
                    self.regPopulation = [self.regPopulation Xmutate(i).mutateordershuffle()];
                elseif r == 2
                    % Memory mutation
                    self.regPopulation = [self.regPopulation Xmutate(i).mutatememory(self.Qnmax,self.Qpmax)];
                elseif r == 3
                    % Order mutation
                    self.regPopulation = [self.regPopulation Xmutate(i).mutateorderincrement()];
                end
                if(self.verbosity>=2) fprintf("%s\n", self.regPopulation(end).print()); end
            end
        end

        %% removerepeated: deletes repeated regs
        function self = removerepeated(self)
            % We add the first three regressors so we ensure they are not
            % lost
            self.regPopulation = [self.regPopulation Regressor([0],[],[])];
            self.regPopulation = [self.regPopulation Regressor([],[0],[])];
            self.regPopulation = [self.regPopulation Regressor([],[],[0])];

            % Canonical form
            for i = 1:length(self.regPopulation)
                self.regPopulation(i).deriveEnvelopeTerms();
                self.regPopulation(i).sortindexes();
            end

            % Removing regressors with order > Pmax
            purgedOrder = 0;
            regdelete = [];
            for i = 1:length(self.regPopulation)
                regOrder = length(self.regPopulation(i).X) + length(self.regPopulation(i).Xconj) + length(self.regPopulation(i).Xenv);
                if regOrder>self.Pmax
                    regdelete = [regdelete i];
                    purgedOrder = purgedOrder + 1;
                end
            end
            self.regPopulation(regdelete) = [];
            if(self.verbosity>=2)
                fprintf('Number of regressors removed (Pmax): %d\n', purgedOrder);
            end

            % Removing regressors with memory > Qmax
            purgedMemory = 0;
            regdelete = [];
            for i = 1:length(self.regPopulation)
                regMemory = max([self.regPopulation(i).X self.regPopulation(i).Xconj self.regPopulation(i).Xenv]);
                if regMemory>self.Qnmax
                    regdelete = [regdelete i];
                    purgedMemory = purgedMemory + 1;
                end
            end
            self.regPopulation(regdelete) = [];
            if(self.verbosity>=2)
                fprintf('Number of regressors removed (Qmax): %d\n', purgedMemory);
            end

            % Removing repeated regressors
            purgedRep = 0;
            for i = 1:length(self.regPopulation)
                j = i+1;
                while j <= length(self.regPopulation)
                    if self.regPopulation(i).equals(self.regPopulation(j))
                        if(self.verbosity>=2)
                            fprintf('Regressor removed: %s\n', self.regPopulation(i).print());
                        end
                        self.regPopulation(j) = [];
                        j = j-1;
                        purgedRep = purgedRep +1;
                    end
                    j = j+1;
                end
            end
            if(self.verbosity>=2)
                fprintf('Number of repeated regressors removed: %d\n', purgedRep);
            end
        end

        %% buildUcustomX
        %   Builds the regressor matrix for a different pair x-y
        function [X,yX] =  buildUcustomX(self, x, y, n)
            nRows = numel(n(1+self.Qpmax:end-self.Qnmax));
            nRegs = length(self.regPopulation);
            estimatedGB = nRows * nRegs * 16 / 2^30;
            fprintf('[regressorManager.buildUcustomX] Materializando U completa: rows=%d regs=%d estimado=%.2f GB\n', ...
                nRows, nRegs, estimatedGB);
            X = complex(zeros(nRows, nRegs));
            for i = 1:length(self.regPopulation)
                self.regPopulation(i).buildRegressor(x, n, self.Qpmax, self.Qnmax);
                X(:, i) = self.regPopulation(i).reg(:);
            end
            yX = y(n(1+self.Qpmax:end-self.Qnmax));
        end

        function [self] =  clearRegressors(self)
            for i = 1:length(self.regPopulation)
                self.regPopulation(i).reg=[];
            end
        end
        %% regress
        %   Regresses the model
        function self = regress(self)
            self.buildX();
            normX = vecnorm(self.X);
            normX(normX == 0) = 1;
            Xn = self.X ./ normX;
            hNorm = Xn \ self.yX;
            ymod = self.X * (hNorm ./ normX(:));
            denom = norm(self.yX,2);
            if denom == 0
                self.nmse = Inf;
            else
                self.nmse = 20*log10(norm(ymod-self.yX,2)/denom);
            end
        end
        %% prepareForSave: Clears unnecessary data to reduce object size before saving
        % This function retains only the essential data in the object,
        % removing the rest to minimize the object's size for efficient storage.
        % It clears the `reg` field in each regressor, as well as the fields
        % `x`, `y`, `n`, and `U`.
        function self = prepareForSave(self)
            % Clear the 'reg' field for each regressor in the population
            for i = 1:length(self.regPopulation)
                self.regPopulation(i).reg = [];
            end

            % Clear unnecessary fields to save space
            self.x = [];
            self.y = [];
            self.n = [];
            self.X = [];
            self.yX = [];
        end
    end
end
