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
            % Initialize the Regressor object's properties
            self.X = X;
            self.Xconj = Xconj;
            self.Xenv = Xenv;
        end
        %% Compares two Regressor objects for equality
        % Returns true if the X, Xconj, and Xenv arrays are equal, or all empty
        function equal = equals(self, other)
            % Check if the arrays are equal or both empty
            equal = (isequal(self.X, other.X) || (isempty(self.X) && isempty(other.X))) && ...
                (isequal(self.Xconj, other.Xconj) || (isempty(self.Xconj) && isempty(other.Xconj))) && ...
                (isequal(self.Xenv, other.Xenv) || (isempty(self.Xenv) && isempty(other.Xenv)));
        end

        %% Generates a compact string representation of the Regressor
        % Formats the X, Xconj, and Xenv terms into a human-readable expression
        function S = print(self)
            S = [];
            self.sortindexes(); % Ensure the arrays are sorted

            % Process the X terms
            [uX] = unique(self.X);
            for ir = 1:length(uX)
                if uX(ir)==0
                    S = sprintf('%s x(n)^%d', S, sum(uX(ir) == self.X));
                else
                    S = sprintf('%s x(n%+d)^%d', S, -uX(ir), sum(uX(ir) == self.X));
                end
            end


            % Process the Xconj terms
            [uXconj] = unique(self.Xconj);
            for ir = 1:length(uXconj)
                if uXconj(ir)==0
                    S = sprintf('%s x*(n)^%d', S, sum(uXconj(ir) == self.Xconj));
                else
                    S = sprintf('%s x*(n%+d)^%d', S, -uXconj(ir), sum(uXconj(ir) == self.Xconj));
                end
            end

            % Process the Xenv terms
            [uXenv] = unique(self.Xenv);
            for ir = 1:length(uXenv)
                if uXenv(ir)==0
                    S = sprintf('%s |x(n)|^%d', S, sum(uXenv(ir) == self.Xenv));
                else
                    S = sprintf('%s |x(n%+d)|^%d', S, -uXenv(ir), sum(uXenv(ir) == self.Xenv));
                end
            end
        end

        %% Sorts the X, Xconj, and Xenv arrays in ascending order
        function self = sortindexes(self)
            if ~isempty(self.X)
                self.X = sort(self.X);
            end
            if ~isempty(self.Xconj)
                self.Xconj = sort(self.Xconj);
            end
            if ~isempty(self.Xenv)
                self.Xenv = sort(self.Xenv);
            end
        end

        %% Derives the Xenv terms from the X and Xconj arrays
        % The relationship X * Xconj = Xenv^2 is used to compute the Xenv terms
        function self = deriveEnvelopeTerms(self)
            % Iterate through the unique elements in the X array
            for i = unique(self.X)
                % Count the occurrences of the current element in X and Xconj
                occurrencesInX = sum(i == self.X);
                occurrencesInXconj = sum(i == self.Xconj);

                % Take the minimum of the occurrences to determine the number of Xenv terms
                numXenvTerms = min([occurrencesInX, occurrencesInXconj]);

                % Add the Xenv terms and remove the corresponding elements from X and Xconj
                if numXenvTerms > 0
                    self.Xenv = [self.Xenv, i * ones(1, 2 * numXenvTerms)];
                    self.X = self.X(setdiff(1:length(self.X), find(i == self.X, numXenvTerms)));
                    self.Xconj = self.Xconj(setdiff(1:length(self.Xconj), find(i == self.Xconj, numXenvTerms)));
                end
            end
        end
    end
end
