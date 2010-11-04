function [data] = processFault(fileBase,numProc,numLevels,tauRef,vr,slope,x0,Z)

ML = [];
for p = 0:(numProc-1)
    for l = 0:(numLevels-1)
        fileName = sprintf('%s.p%d.l%d',fileBase,p,l);
        tmp = load(fileName);
        ML = [ML;tmp];
    end
end

data.time      = ML(1,1);
data.numLevels = numLevels;
data.tauRef = tauRef;
data.vr = vr;
data.pulseEdge = data.time * data.vr;
data.slope = slope;
L = ML(:,2);
for lvl = 1:numLevels
    I     = find(L == (lvl-1));
    ind   = ML(I, 3);
    [x,J] = sort(ML(I, 4));
    I     = I(J);
    vx    = ML(I, 5);
    vy    = ML(I, 6);
    vz    = ML(I, 7);
    sxx   = ML(I, 8);
    syy   = ML(I, 9);
    szz   = ML(I,10);
    sxy   = ML(I,11);
    sxz   = ML(I,12);
    syz   = ML(I,13);

    wpx  = -sxy - Z*vx;
    wpz  = -syz - Z*vz;
    dx = 100/(128*4^(lvl-1));
    tau0 = tauRef +slope*max(0,abs(x+dx/2-x0) - data.pulseEdge);
    tau  = sqrt((wpx+wpx).^2 + (wpz+wpz).^2)/2;
    beta = min(1,2*tau0./tau - 1);
    wmx  = beta.*wpx;
    wmz  = beta.*wpz;
    tau  = sqrt((wmx+wpx).^2 + (wmz+wpz).^2)/2;
    V    = sqrt((wmx-wpx).^2 + (wmz-wpz).^2)/(2*Z);


    data.levels(lvl).ind = ind;
    data.levels(lvl).x   = x  ;
    data.levels(lvl).vx  = vx ;
    data.levels(lvl).vy  = vy ;
    data.levels(lvl).vz  = vz ;
    data.levels(lvl).sxx = sxx;
    data.levels(lvl).syy = syy;
    data.levels(lvl).szz = szz;
    data.levels(lvl).sxy = sxy;
    data.levels(lvl).sxz = sxz;
    data.levels(lvl).syz = syz;
    data.levels(lvl).tau = tau;
    data.levels(lvl).V   = V;
    data.levels(lvl).tau0 = tau0;
end

x   = [];
vx  = [];
vy  = [];
vz  = [];
sxx = [];
syy = [];
szz = [];
sxy = [];
sxz = [];
syz = [];
tau = [];
V   = [];
L   = [];
for lvl = numLevels:-1:1
    [~,I] = setdiff(data.levels(lvl).x,x);
    x   = [x  ;data.levels(lvl).x(I)  ];
    L   = [L  ; I*0 + lvl];
    vx  = [vx ;data.levels(lvl).vx(I) ];
    vy  = [vy ;data.levels(lvl).vy(I) ];
    vz  = [vz ;data.levels(lvl).vz(I) ];
    sxx = [sxx;data.levels(lvl).sxx(I)];
    syy = [syy;data.levels(lvl).syy(I)];
    szz = [szz;data.levels(lvl).szz(I)];
    sxy = [sxy;data.levels(lvl).sxy(I)];
    sxz = [sxz;data.levels(lvl).sxz(I)];
    syz = [syz;data.levels(lvl).syz(I)];
    tau = [tau;data.levels(lvl).tau(I)];
    V   = [V  ;data.levels(lvl).V(I)  ];
end

[~,I] = sort(x);
data.x   = x(I); 
data.vx  = vx(I); 
data.vy  = vy(I); 
data.vz  = vz(I); 
data.sxx = sxx(I);
data.syy = syy(I);
data.szz = szz(I);
data.sxy = sxy(I);
data.sxz = sxz(I);
data.syz = syz(I);
data.tau = tau(I);
data.V   = V(I);
data.L   = L(I);
