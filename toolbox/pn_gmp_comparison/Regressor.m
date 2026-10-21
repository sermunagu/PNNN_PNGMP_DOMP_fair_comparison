classdef Regressor < handle
    properties
        X        % Memory taps (lags) for the original input signal.
        Xconj    % Memory taps (lags) for the complex conjugate of the input signal.
        Xenv     % Memory taps (lags) for the absolute value (envelope) of the input signal.
    end
    %%
    methods
        % Builder, Inputs:
        %   X = [m0,m1,m2,m3...]
        %   Xconj = [m0,m1,m2,m3...]
        %   Xenv = [m0,m1,m2,m3...]
        function self = Regressor(X, Xconj, Xenv)
            self.X = X;
            self.Xconj = Xconj;
            self.Xenv = Xenv;
        end
        %% Compares two Regressor objects for equality
        % Returns true if the X, Xconj, and Xenv arrays are equal, or all empty
        function equal = equals(self, other)
            equal = (isequal(self.X, other.X) || (isempty(self.X) && isempty(other.X))) && ...
                (isequal(self.Xconj, other.Xconj) || (isempty(self.Xconj) && isempty(other.Xconj))) && ...
                (isequal(self.Xenv, other.Xenv) || (isempty(self.Xenv) && isempty(other.Xenv)));
        end

        %% Sorts the X, Xconj, and Xenv arrays in ascending order
        function self = sortLags(self)
            self.X = sort(self.X);
            self.Xconj = sort(self.Xconj);
            self.Xenv = sort(self.Xenv);
        end

        %% Derives the Xenv terms from the X and Xconj arrays
        % The relationship X * Xconj = Xenv^2 is used to compute the Xenv terms
        function self = deriveEnvelopeTerms(self)
            for i = unique(self.X)
                occurrencesInX = sum(i == self.X);
                occurrencesInXconj = sum(i == self.Xconj);
                numXenvTerms = min([occurrencesInX, occurrencesInXconj]);
                if numXenvTerms > 0
                    self.Xenv = [self.Xenv, i * ones(1, 2 * numXenvTerms)];
                    self.X = self.X(setdiff(1:length(self.X), find(i == self.X, numXenvTerms)));
                    self.Xconj = self.Xconj(setdiff(1:length(self.Xconj), find(i == self.Xconj, numXenvTerms)));
                end
            end
        end
    end
end
