function resolved = resolveLatexmkCommand(preferredCommand)
% resolveLatexmkCommand - Resolve latexmk, including MATLAB-Perl on MiKTeX.
% MiKTeX may install latexmk.pl without registering a Perl script engine.
% In that case the Perl already shipped with MATLAB is used without changing
% PATH, the MiKTeX configuration, or any system installation.

if nargin < 1 || isempty(preferredCommand)
    preferredCommand = 'latexmk';
end
nativePrefix = quoteArgument(char(string(preferredCommand)));
[status, versionText] = system(char(nativePrefix + " --version"));
if status == 0
    resolved = struct('commandPrefix', string(nativePrefix), ...
        'versionText', string(versionText), 'usedMatlabPerlFallback', false);
    return;
end

if ispc
    perlExecutable = fullfile(matlabroot, 'sys', 'perl', 'win32', ...
        'bin', 'perl.exe');
    latexmkScript = fullfile(getenv('LOCALAPPDATA'), 'Programs', ...
        'MiKTeX', 'scripts', 'latexmk', 'latexmk.pl');
    if isfile(perlExecutable) && isfile(latexmkScript)
        fallbackPrefix = quoteArgument(perlExecutable) + " " + ...
            quoteArgument(latexmkScript);
        [fallbackStatus, fallbackText] = system( ...
            char(fallbackPrefix + " --version"));
        if fallbackStatus == 0
            resolved = struct('commandPrefix', fallbackPrefix, ...
                'versionText', string(fallbackText), ...
                'usedMatlabPerlFallback', true);
            return;
        end
        versionText = string(versionText) + newline + string(fallbackText);
    end
end

error('resolveLatexmkCommand:Unavailable', ...
    'latexmk could not be executed: %s', strtrim(versionText));
end

function value = quoteArgument(value)
quote = string(char(34));
value = quote + replace(string(value), quote, quote + quote) + quote;
end
