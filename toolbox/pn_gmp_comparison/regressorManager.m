classdef regressorManager < matlab.mixin.Copyable
    properties
        regPopulation  % Regressors vector
        Qpmax          % Maximum noncausal memory.
        Qnmax          % Maximum causal memory.
        Pmax           % Maximum nonlinear order.
        % Output config:
        verbosity      % Level of output detail.
    end
    %%
    methods
        %% regressorManager: builder
        function self = regressorManager(~,~,config)
            self.Qpmax = config.Qpmax;
            self.Qnmax = config.Qnmax;
            self.Pmax = config.Pmax;
            self.regPopulation = [];
            self.verbosity = 0;
        end
        %% initialization: creates initial population
        %       Create the population with a defined initial set:
        function self = initialization(self)
            modelsConfig;

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

    end
end
