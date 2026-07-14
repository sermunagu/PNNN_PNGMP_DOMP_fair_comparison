classdef Regressor < handle
    properties
        X        % Memory taps (lags) for the original input signal.
        Xconj    % Memory taps (lags) for the complex conjugate of the input signal.
        Xenv     % Memory taps (lags) for the absolute value (envelope) of the input signal.
        score    % Fitness score associated with this regressor.
        reg      % Regressor matrix constructed using X, Xconj, and Xenv.
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
            self.score = 0;
            self.reg = [];
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

        %% mutatememory: Creates a mutated copy of the current regressor
        % This function generates a mutated version of the current regressor by
        % randomly selecting and incrementing one memory tap in the arrays `X`, `Xconj`,
        % or `Xenv`, which represent different components of the regressor.
        %
        function newreg = mutatememory(self,Qnmax,Qpmax)
            newreg = Regressor(self.X,self.Xconj,self.Xenv);
            ind = randi([1 length(self.X)+length(self.Xconj)+length(self.Xenv)]);
			
			memorysign = +1;
			Pneg = 0.1; % Probability of reducing one memory tap 
			if rand<Pneg
				memorysign = -1;
			end
			
            if ind <= length(self.X)
                newreg.X(ind) = max(min(newreg.X(ind) + memorysign,Qpmax),-Qnmax);
            elseif ind <= length(self.X)+length(self.Xconj)
                ind = ind-length(self.X);
                newreg.Xconj(ind) = max(min(newreg.Xconj(ind) + memorysign,Qpmax),-Qnmax);
            else
                ind = ind - length(self.X) - length(self.Xconj);
                newreg.Xenv(ind) = max(min(newreg.Xenv(ind) + memorysign,Qpmax),-Qnmax);
            end
        end

        %% mutateordershuffle: Alters the order of terms in the regressor by moving elements between X, Xconj, and Xenv
        % This function generates a mutated version of the current regressor by
        % randomly moving one element between the arrays `X`, `Xconj`, and `Xenv`,
        % which represent different components of the regressor.
        %
        % The mutation process proceeds as follows:
        %   - A new `Regressor` object `newreg` is created as a copy of the current
        %     regressor, initialized with the values of `X`, `Xconj`, and `Xenv`.
        %   - A random index `ind` is selected from the combined length of `X`,
        %     `Xconj`, and `Xenv`, determining which array will have an element moved.
        %   - Based on the range of `ind`, an element is randomly selected from one
        %     of the arrays (`X`, `Xconj`, or `Xenv`):
        %       - **If `ind` falls within the range of `X`**:
        %           - Randomly select an element from `X`.
        %           - Move this element to either `Xconj` or `Xenv` with a 50% chance
        %             for each option.
        %           - Remove the selected element from `X`.
        %       - **If `ind` falls within the range of `Xconj`**:
        %           - Randomly select an element from `Xconj`.
        %           - Move this element to either `X` or `Xenv` with a 50% chance
        %             for each option.
        %           - Remove the selected element from `Xconj`.
        %       - **If `ind` falls within the range of `Xenv`**:
        %           - Randomly select an element from `Xenv`.
        %           - Move this element to either `X` or `Xconj` with a 50% chance
        %             for each option.
        %           - Remove the selected element from `Xenv`.
        %
        % Output:
        %   - newreg: The mutated regressor, with one element moved between `X`,
        %     `Xconj`, or `Xenv`.
        %
        % This mutation introduces diversity in the regressor population by reordering
        % elements between the different arrays, which can help the algorithm explore
        % new configurations of terms for better solution discovery.
        function newreg = mutateordershuffle(self)
            newreg = Regressor(self.X,self.Xconj,self.Xenv);
            ind = randi([1 length(self.X)+length(self.Xconj)+length(self.Xenv)]);
            if ind <= length(self.X)
                selected = randi([1,length(self.X)]);
                if rand<0.5
                    newreg.Xconj = [newreg.Xconj,self.X(selected)];
                else
                    newreg.Xenv = [newreg.Xenv,self.X(selected)];
                end
                newreg.X(selected) = [];
            elseif ind <= length(self.Xconj)+length(self.X)
                selected = randi([1,length(self.Xconj)]);
                if rand<0.5
                    newreg.X = [newreg.X,self.Xconj(selected)];
                else
                    newreg.Xenv = [newreg.Xenv,self.Xconj(selected)];
                end
                newreg.Xconj(selected) = [];
            else
                selected = randi([1,length(self.Xenv)]);
                if rand<0.5
                    newreg.X = [newreg.X,self.Xenv(selected)];
                else
                    newreg.Xconj = [newreg.Xconj,self.Xenv(selected)];
                end
                newreg.Xenv(selected) = [];
            end
        end
        %% mutateorderincrement: Adds a new zero term to the regressor's X or Xconj arrays
        % This function generates a mutated version of the current regressor by
        % appending a zero to either the `X` or `Xconj` arrays with a 50% chance
        % for each option. This mutation introduces an additional term into the
        % regressor, potentially altering its structure.
        %
        % The mutation process proceeds as follows:
        %   - A new `Regressor` object `newreg` is created as a copy of the current
        %     regressor, initialized with the current values of `X`, `Xconj`, and `Xenv`.
        %   - A random decision is made (using `rand < 0.5`) to either:
        %       - Append a zero to the `X` array if the random value is less than 0.5.
        %       - Append a zero to the `Xconj` array if the random value is 0.5 or greater.
        %
        % Output:
        %   - newreg: The mutated regressor, which includes an additional zero term in
        %     either `X` or `Xconj`.
        %
        % This mutation adds flexibility to the regressor population by allowing
        % growth in the regressor's representation. It provides the potential for new
        % configurations, especially when zeros represent a placeholder or initial
        % value that can later be modified or replaced in future mutations.
        function newreg = mutateorderincrement(self)
            newreg = Regressor(self.X,self.Xconj,self.Xenv);
            if rand<0.5
                newreg.X = [newreg.X 0];
            else
                newreg.Xconj = [newreg.Xconj 0];
            end
        end

        %% crossover: Performs crossover between two parents to create a new regressor
        % This function implements the crossover operation in a genetic algorithm
        % to combine the genetic material of two parent regressors and create a
        % new offspring. The crossover is performed separately on the three
        % components of the regressor: `X`, `Xconj`, and `Xenv`.
        function newreg = crossover(self, parent)
            Xn = [self.X parent.X];
            Xconjn = [self.Xconj parent.Xconj];
            Xenvn = [self.Xenv parent.Xenv];
            newreg = Regressor(Xn,Xconjn,Xenvn);
        end
        function newreg = crossoverkk(self, parent)
            ph = ([self.X parent.X]);
            ph = ph(randperm(length(ph)));
            Xn = ph(1:floor((length(self.X)+length(parent.X))/2));
            ph = ([self.Xconj parent.Xconj]);
            ph = ph(randperm(length(ph)));
            Xconjn = ph(1:floor((length(self.Xconj)+length(parent.Xconj))/2));
            ph = ([self.Xenv parent.Xenv]);
            ph = ph(randperm(length(ph)));
            Xenvn = ph(1:floor((length(self.Xenv)+length(parent.Xenv))/2));
            newreg = Regressor(Xn,Xconjn,Xenvn);
        end
        %% buildRegressor: Creates the regressor matrix for model building
        % This function constructs the regressor matrix `reg` used in a predictive or
        % estimation model by combining multiple transformations of the input signal `x`.
        % The matrix includes:
        %   - Original values of `x` (indexed by the vector `X`)
        %   - Complex conjugates of `x` (indexed by the vector `Xconj`)
        %   - Absolute values of `x` (indexed by the vector `Xenv`)
        %
        % Parameters:
        %   self    - The object instance of the class (modified to include the regressor)
        %   x       - Input signal vector
        %   n       - Index vector for sampling positions
        %   Qpmax   - Maximum positive lag (determines start point for indexing)
        %   Qnmax   - Maximum negative lag (determines end point for indexing)
        %
        % The regressor `self.reg` is initialized to a vector of ones, then updated by
        % element-wise multiplication with different transformations of `x`, following
        % these steps:
        %   1. Multiply by the values of `x` specified in `self.X`.
        %   2. Multiply by the complex conjugate values of `x` specified in `self.Xconj`.
        %   3. Multiply by the absolute values of `x` specified in `self.Xenv`.
        %
        % Output:
        %   The function updates `self.reg` to store the final regressor matrix, with each
        %   element reflecting the combined transformations of `x` as specified by `X`,
        %   `Xconj`, and `Xenv`.
        function self = buildRegressor(self, x, n, Qpmax, Qnmax)
            reg = ones(size(n(1+Qpmax:end-Qnmax)));
            for i=1:length(self.X)
                reg = reg.*x(n(1+Qpmax-self.X(i):end-Qnmax-self.X(i)));
            end
            for i=1:length(self.Xconj)
                reg = reg.*conj(x(n(1+Qpmax-self.Xconj(i):end-Qnmax-self.Xconj(i))));
            end
            for i=1:length(self.Xenv)
                reg = reg.*abs(x(n(1+Qpmax-self.Xenv(i):end-Qnmax-self.Xenv(i))));
            end
            self.reg = reg;
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
