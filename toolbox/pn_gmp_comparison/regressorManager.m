classdef regressorManager < matlab.mixin.Copyable
    properties
        regPopulation  % Regressors vector
        Qpmax          % Maximum noncausal memory.
        Qnmax          % Maximum causal memory.
        Pmax           % Maximum nonlinear order.
    end
    %%
    methods
        %% regressorManager: builder
        function self = regressorManager(~,~,config)
            self.Qpmax = config.Qpmax;
            self.Qnmax = config.Qnmax;
            self.Pmax = config.Pmax;
            self.regPopulation = [];
        end
        %% buildPopulation: creates initial population
        %       Create the population with a defined initial set:
        function self = buildPopulation(self)
            % Deterministic GMP population dimensions.
            Pgmp = 13; Lgmp = 10; Mgmp = 2;

            PgmpEff = max(1, min(Pgmp, self.Pmax));
            KaEff = 0:(PgmpEff-1);
            KbcEff = 1:(PgmpEff-1);

            for k = 1:length(KaEff)
                for l = 0:min(Lgmp, self.Qpmax)
                    X = l;
                    Xenv = repmat(l,1,KaEff(k));

                    self.regPopulation = [self.regPopulation Regressor(X,[],Xenv)];
                end
            end

            for k = 1:length(KbcEff)
                for l = 0:min(Lgmp, self.Qpmax)
                    mMax = min(Mgmp, self.Qpmax - l);
                    for m = 1:mMax

                        X = l;
                        Xenv = repmat(l+m,1,KbcEff(k));

                        self.regPopulation = [self.regPopulation Regressor(X,[],Xenv)];
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
                    end
                end
            end
        end
        %% Add auxiliaries, canonicalize, filter limits and remove duplicate regressors
        function self = finalizePopulation(self)
            % We add the first three regressors so we ensure they are not
            % lost
            self.regPopulation = [self.regPopulation Regressor(0,[],[])];
            self.regPopulation = [self.regPopulation Regressor([],0,[])];
            self.regPopulation = [self.regPopulation Regressor([],[],0)];

            % Canonical form
            for i = 1:length(self.regPopulation)
                self.regPopulation(i).deriveEnvelopeTerms();
                self.regPopulation(i).sortLags();
            end

            % Removing regressors with order > Pmax
            regdelete = false(1, numel(self.regPopulation));
            for i = 1:length(self.regPopulation)
                regOrder = length(self.regPopulation(i).X) + length(self.regPopulation(i).Xconj) + length(self.regPopulation(i).Xenv);
                regdelete(i) = regOrder > self.Pmax;
            end
            self.regPopulation(regdelete) = [];

            % Removing regressors with memory > Qmax
            regdelete = false(1, numel(self.regPopulation));
            for i = 1:length(self.regPopulation)
                regMemory = max([self.regPopulation(i).X self.regPopulation(i).Xconj self.regPopulation(i).Xenv]);
                regdelete(i) = regMemory > self.Qnmax;
            end
            self.regPopulation(regdelete) = [];

            % Removing repeated regressors
            for i = 1:length(self.regPopulation)
                j = i+1;
                while j <= length(self.regPopulation)
                    if self.regPopulation(i).equals(self.regPopulation(j))
                        self.regPopulation(j) = [];
                        j = j-1;
                    end
                    j = j+1;
                end
            end
        end

    end
end
