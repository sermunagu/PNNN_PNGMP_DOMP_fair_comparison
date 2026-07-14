% GVG configuration for initialization

% Full Volterra
Pfv=13;
Mfv=3;

% CVS
Pcvs=13;
Mcvs=3;

% MP
Pmp = 13;                
Mmp = 10;

% DDR
Pddr = 13;                
Mddr = 10;

%GMP
Pgmp = 13;                
Lgmp = 10;
Mgmp = 2;
Ka = [0:(Pgmp-1)];
La = Lgmp*ones(size(Ka));
Kb = [1:(Pgmp-1)];
Lb = Lgmp*ones(size(Kb));
Mb = Mgmp*ones(size(Kb));
Kc = [1:(Pgmp-1)];
Lc = Lgmp*ones(size(Kc));
Mc = Mgmp*ones(size(Kc));