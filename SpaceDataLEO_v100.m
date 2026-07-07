%% LEO_Constellation_ISL_Simulation_4terminals.m
%
% Time-evolution simulation of a Walker-Delta LEO satellite constellation
% with:
%   P  orbital planes, equally spaced in RAAN
%   Ns satellites per plane (orbital ring), equally spaced in mean anomaly
%   inclination theta
%
% TERMINAL MODEL: each satellite has exactly FOUR ISL terminals:
%   Terminal 1 (PERMANENT)    -> previous satellite in the same ring
%   Terminal 2 (PERMANENT)    -> next satellite in the same ring
%   Terminal 3 (INTERMITTENT) -> one confirmed bidirectional link toward pPrev
%   Terminal 4 (INTERMITTENT) -> one confirmed bidirectional link toward pNext
%
% INTER-RING LINK ESTABLISHMENT (three-pass procedure each timestep):
%
%   Pass 1 - proposals (range + hemispherical FoR only):
%     Every terminal 3 and terminal 4 independently evaluates ALL Ns
%     satellites in its designated adjacent plane and proposes a link to
%     every one that satisfies BOTH:
%       (a) hemispherical FoR: dot(r2-r1, boresight) > 0
%       (b) range:             norm(r2-r1) <= maxISLrange
%     Earth occlusion and Sun-blinding are NOT applied here.
%     The result is a proposal SET (not a single choice) per terminal:
%       propSet3{p,s} - satellites in pPrev proposed by terminal-3
%       propSet4{p,s} - satellites in pNext proposed by terminal-4
%
%   Pass 2 - mutual-pair identification:
%     A pair ((p,s) via terminal-4, (pNext,s2) via terminal-3) is MUTUAL if:
%       s2 in propSet4{p,s}   AND   s in propSet3{pNext,s2}
%     All mutual pairs are collected with their Euclidean distance.
%
%   Pass 3 - greedy nearest-first assignment with Earth + Sun gates:
%     Mutual pairs are sorted by ascending distance. Each pair is accepted
%     (link confirmed) if and only if:
%       (i)  terminal-4 of (p,s)    is still free, AND
%       (ii) terminal-3 of (pNext,s2) is still free, AND
%       (iii) the link is not Earth-occluded, AND
%       (iv)  the link is not Sun-blinded.
%     Once a terminal is assigned to one confirmed link it is marked busy
%     and cannot accept another link in the same timestep (at most one
%     confirmed link per terminal, as required by the 4-terminal model).
%
% Each confirmed link contributes 2*T_I to aggregate capacity (T_I per
% endpoint terminal) and credits T_I to each endpoint plane in the
% bisection calculation.
%
% Orbits are propagated as simple two-body circular Keplerian orbits
% (good approximation for near-circular LEO constellations such as
% Starlink, e~0).
%
% Author: Alberto Perotti @ ASTRACOM/Politecnico di Torino
% (alberto.perotti@polito.it)
% (based on Claude output)

clear; close all; clc;

%% ---------------- USER PARAMETERS ----------------
P        = 20;          % number of orbital planes
Ns       = 10;          % satellites per plane (per ring)
theta    = 53;          % inclination [deg]
altitude = 600;         % altitude [km] (circular orbit)
F        = 10;          % Walker-Delta phasing factor (integer, 0..P-1)

T_P = 100;      % permanent (intra-plane, ring) terminal throughput [Gbps]
T_I = 100;      % intermittent (inter-plane) terminal throughput [Gbps]

maxISLrange = 3000;    % [km] max range for a feasible cross-plane ISL

% FIELD OF REGARD (hemispherical, opposite boresights for terminals 3 & 4)
% Each intermittent terminal has a hemispherical FoR whose boresight axis
% is the orbit-normal of the hosting satellite:
%   n_p = Rz(Om_p) * Rx(i) * z_hat
%
%   Terminal 3 boresight: +n_p   Terminal 4 boresight: -n_p
%
% KEY GEOMETRY: the sign of dot(r_q - r_p, n_p) changes over the orbit.
%   dot(dv, n_p) ≈ -a·DeltaOmega·cos(u)·sin(i)
%   u in [0°, 90°) -> dot < 0 -> pNext in -n_p hemi (T4 side)
%   u in (90°,270°) -> dot > 0 -> pNext in +n_p hemi (T3 side)
% A fixed plane designation (T3->pPrev only, T4->pNext only) would leave
% both terminals idle for half the orbit (the half when the target plane
% is on the wrong side of the boresight). FIX: each terminal searches
% BOTH adjacent planes; the FoR naturally selects whichever plane is
% currently in its hemisphere, giving valid links throughout the orbit.

sunLon_deg = 35;         % Sun longitude: right-ascension of the Sun measured
                          % in the Earth equatorial plane from the X-axis [deg]
sunLat_deg =  0;         % Sun latitude: declination of the Sun above/below the
                          % equatorial plane [deg].  Range -23.4° to +23.4° for
                          % realistic seasonal values; 0° = equinox.
phi_deg    = 30;          % [deg] angular exclusion half-cone around the
                          % Earth-Sun line: an ISL is unavailable if its
                          % line-of-sight forms an angle smaller than phi
                          % with the Earth-Sun line AND both of its
                          % endpoint satellites are currently in sunlight
                          % (not eclipsed by Earth's shadow). A satellite
                          % in Earth's shadow has no direct Sun in its
                          % field of view, so its terminals cannot be
                          % blinded by it.



%% ---------------- PHYSICAL CONSTANTS --------------
Re = 6378.137;             % Earth mean radius [km]
mu = 398600.4418;          % Earth gravitational parameter [km^3/s^2]
a  = Re + altitude;              % semi-major axis (circular orbit) [km]
incl = deg2rad(theta);

Torb   = 2*pi*sqrt(a^3/mu); % orbital period [s]
n_mean = 2*pi/Torb;         % mean motion [rad/s]
c_light = 299792.458;       % speed of light [km/s]

dt      = 10;               % simulation time step [s]
simTime = ceil(Torb/dt) * dt;   % total simulated time [s]

fprintf('Semi-major axis : %.1f km\n', a);
fprintf('Orbital period  : %.2f min\n', Torb/60);

%% ---------------- SUN GEOMETRY ---------------------
% Sun assumed fixed (simulation time << 1 year) at zero declination
% (equatorial plane), arbitrary right-ascension sunLon_deg.
% Sun direction unit vector in ECI, combining longitude and latitude:
%   sunDir = [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)]
sunDir = [cosd(sunLat_deg)*cosd(sunLon_deg);
          cosd(sunLat_deg)*sind(sunLon_deg);
          sind(sunLat_deg)];             % unit vector Earth->Sun
phi    = deg2rad(phi_deg);

%% ---------------- WALKER-DELTA GEOMETRY -----------
dRAAN  = 2*pi/P;                       % RAAN spacing between planes
RAAN0  = (0:P-1)*dRAAN;                % [rad]
dPhase = 2*pi*F/(P*Ns);                % inter-plane phase offset

M0 = zeros(P,Ns);                      % initial mean anomaly [rad]
for p = 1:P
    for s = 1:Ns
        M0(p,s) = 2*pi*(s-1)/Ns + (p-1)*dPhase;
    end
end

%% ---------------- TIME LOOP SETUP -----------------
tVec = 0:dt:simTime;
Nt   = length(tVec);

capacityPerm   = zeros(Nt,1);
capacityInterm = zeros(Nt,1);
capacityTotal  = zeros(Nt,1);
numIntermLinks    = zeros(Nt,1);  % # confirmed bidirectional inter-ring links
numTotalProposals = zeros(Nt,1);  % # proposals made (range + FoR gate only)
numMutualPairs    = zeros(Nt,1);  % # mutual pairs found before Earth/Sun gates
numHemiRejected   = zeros(Nt,1);  % # candidates failing hemispherical FoR
numEarthBlocked   = zeros(Nt,1);  % # mutual pairs rejected by Earth occlusion
numSunBlockedInterm = zeros(Nt,1);% # mutual pairs rejected by Sun-blinding
numBusyTerminal   = zeros(Nt,1);  % # mutual pairs skipped: terminal already used
numSunBlockedPerm = zeros(Nt,1);  % # permanent links blinded by the Sun
numEclipsed       = zeros(Nt,1);  % # satellites currently in Earth's shadow
bisectionBW  = zeros(Nt,1);
oversubRatio = zeros(Nt,1);

% --- Latency metrics ---
% meanHopLatency : mean single-hop propagation delay across ALL active links [ms]
% meanE2ELatency : mean end-to-end shortest-path delay across reachable pairs [ms]
% maxE2ELatency  : network diameter = max shortest-path delay [ms]
meanHopLatency = zeros(Nt,1);
meanE2ELatency = zeros(Nt,1);
maxE2ELatency  = zeros(Nt,1);

% --- Reliability metrics ---
% netConnected   : true if all satellite pairs are mutually reachable
% fracLinksActive: (active links) / (max possible links = 2*P*Ns)
% fracReachable  : fraction of ordered satellite pairs (i!=j) with finite path
netConnected    = false(Nt,1);
fracLinksActive = zeros(Nt,1);
fracReachable   = zeros(Nt,1);

% --- Hop-count metrics (unweighted shortest paths) ---
% meanHops : mean hop count across all reachable satellite pairs
% maxHops  : maximum hop count (diameter in hops) across reachable pairs
meanHops = zeros(Nt,1);
maxHops  = zeros(Nt,1);

% --- Solar blinding step counters per terminal ---
% N_nodes defined here so counter arrays can be pre-allocated before the loop.
N_nodes = P*Ns;
% permBlindCount(kk)   : number of timesteps ring link kk is sun-blinded.
%   Ring link kk = permLinksAll(kk,:) = [p, s1, s2] corresponds to
%   terminal-2 of satellite s1 and terminal-1 of satellite s2 in plane p;
%   both terminals blind/unblind together, so one counter per link captures
%   the blinding exposure of the associated terminal pair.
% termBlindCount3(p,s) : number of timesteps terminal-3 of satellite (p,s)
%   is blocked by Sun-blinding in Pass 3 (its nearest mutual-pair candidate
%   is rejected because the link falls within phi of the Earth-Sun line
%   while both endpoints are in sunlight).
% termBlindCount4(p,s) : same for terminal-4.
permBlindCount  = zeros(N_nodes, 1);  % intra-ring,  one entry per ring link
termBlindCount3 = zeros(P, Ns);       % inter-ring terminal-3, P x Ns
termBlindCount4 = zeros(P, Ns);       % inter-ring terminal-4, P x Ns

% --- Solar blinding step counters (NEW) ---
% N_nodes defined here so the counter arrays can be pre-allocated.
N_nodes = P*Ns;
% permBlindCount(kk)   : number of timesteps ring link kk is sun-blinded.
%   Each ring link kk in permLinksAll = [p, s1, s2] corresponds to
%   terminal-2 of satellite s1 and terminal-1 of satellite s2 in plane p;
%   both terminals are blinded/available together, so one counter per link
%   captures the blinding exposure of the associated terminal pair.
% termBlindCount3(p,s) : timesteps terminal-3 of satellite (p,s) is blocked
%   by Sun-blinding in Pass 3 (i.e., its best mutual-pair candidate is
%   rejected because the link angle is within phi of the Earth-Sun line
%   while both endpoints are in sunlight).
% termBlindCount4(p,s) : same for terminal-4.
permBlindCount  = zeros(N_nodes, 1);  % intra-ring,  one entry per ring link
termBlindCount3 = zeros(P, Ns);       % inter-ring terminal-3, P×Ns
termBlindCount4 = zeros(P, Ns);       % inter-ring terminal-4, P×Ns

posAll = zeros(P,Ns,3);   % ECI positions, reused every step
planeNormal = zeros(P,3); % orbit-normal (cross-track) unit vector per plane,
                           % reused every step (constant given Om, incl)

%% ---------------- BISECTION CUT GEOMETRY ----------
% The plane index forms a ring (each plane's intermittent terminals only
% reach its immediate previous/next plane). A bisection that splits the
% planes into two CONTIGUOUS halves (P/2 planes each) therefore only cuts
% intermittent (terminal 3/4) links -- permanent (intra-plane) links
% never cross such a cut. There are P possible rotations of where to
% place the boundary; for each time step we scan all of them and keep
% the one giving the MINIMUM cut capacity, i.e. the true (worst-case)
% bisection bandwidth.
if mod(P,2) ~= 0
    error('P must be even to form two equal-size halves for the bisection cut.');
end
halfP = P/2;
maskGroupA = false(P,P);     % maskGroupA(r,p) = true if plane p is in
                              % "Group A" for rotation r (r = 1..P)
for r = 1:P
    for p = 1:P
        maskGroupA(r,p) = mod(p-1-(r-1), P) < halfP;
    end
end

%% ---------------- FIGURE SETUP --------------------
figure('Color','w','Position',[80 80 950 800]);
hAx = axes; hold(hAx,'on'); axis(hAx,'equal'); grid(hAx,'on');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(45,25);
Lmax = 1.3*a;
xlim([-Lmax Lmax]); ylim([-Lmax Lmax]); zlim([-Lmax Lmax]);

% Earth sphere
[xe,ye,ze] = sphere(40);
surf(hAx, Re*xe, Re*ye, Re*ze, 'FaceColor',[0.3 0.5 0.9], ...
     'EdgeColor','none', 'FaceAlpha',0.35);

hSat  = plot3(hAx, NaN,NaN,NaN, 'k.', 'MarkerSize',10);
hSatEclipsed = plot3(hAx, NaN,NaN,NaN, '.', 'Color',[0.5 0.5 0.5], 'MarkerSize',10);
hPerm = gobjects(0);
hInt  = gobjects(0);
hPermBlocked = gobjects(0);

% Sun direction indicator (arrow pointing outward from Earth, far field)
sunArrowLen = 1.15*Lmax;
quiver3(hAx, 0,0,0, sunDir(1)*sunArrowLen, sunDir(2)*sunArrowLen, sunDir(3)*sunArrowLen, ...
        0, 'Color',[1 0.6 0],'LineWidth',2.5,'MaxHeadSize',0.5);
text(hAx, sunDir(1)*sunArrowLen, sunDir(2)*sunArrowLen, sunDir(3)*sunArrowLen, ...
     '  Sun', 'Color',[1 0.6 0], 'FontWeight','bold');

%% ---------------- MAIN TIME-EVOLUTION LOOP --------
for it = 1:Nt
    t = tVec(it);

    % --- propagate every satellite (circular two-body orbit) ---
    for p = 1:P
        Om = RAAN0(p);
        Rz_Om = [cos(Om) -sin(Om) 0; sin(Om) cos(Om) 0; 0 0 1];
        Rx_i  = [1 0 0; 0 cos(incl) -sin(incl); 0 sin(incl) cos(incl)];
        planeNormal(p,:) = (Rz_Om*Rx_i*[0;0;1])';  % orbit-normal (cross-track)
                                                     % unit vector, constant
                                                     % per plane (circular orbit)
        for s = 1:Ns
            u = M0(p,s) + n_mean*t;        % argument of latitude
            r_pf  = a*[cos(u); sin(u); 0]; % position in orbital plane
            posAll(p,s,:) = Rz_Om*Rx_i*r_pf;
        end
    end

    % --- eclipse status of every satellite (for plotting / stats) ------
    litMask = false(P,Ns);
    for p = 1:P
        for s = 1:Ns
            litMask(p,s) = inSunlight(squeeze(posAll(p,s,:)), sunDir, Re);
        end
    end
    numEclipsed(it) = sum(~litMask(:));

    % --- PERMANENT TERMINALS (1 & 2): ring topology within each plane --
    % Each satellite's terminal 1 connects to the previous satellite in
    % its ring; terminal 2 to the next. Each is tested independently
    % against the Sun-exclusion cone; while blinded, that terminal's
    % throughput drops to zero (link remains topologically permanent).
    permLinksAll = zeros(P*Ns,3);  % [plane, sat(terminal2 owner), nextSat]
    k = 0;
    for p = 1:P
        for s = 1:Ns
            s2 = mod(s,Ns)+1;           % next satellite, wraps around
            k = k+1;
            permLinksAll(k,:) = [p s s2];   % represents BOTH s's terminal-2
                                              % and s2's terminal-1 (same
                                              % physical fiber/laser link)
        end
    end
    permSunBlocked = false(size(permLinksAll,1),1);
    for kk = 1:size(permLinksAll,1)
        p = permLinksAll(kk,1); s1 = permLinksAll(kk,2); s2 = permLinksAll(kk,3);
        r1 = squeeze(posAll(p,s1,:)); r2 = squeeze(posAll(p,s2,:));
        permSunBlocked(kk) = sunBlinded(r1, r2, sunDir, phi, Re);
    end
    permLinks = permLinksAll(~permSunBlocked,:);   % only AVAILABLE perm links
    numSunBlockedPerm(it) = sum(permSunBlocked);
    permBlindCount = permBlindCount + double(permSunBlocked);  % accumulate per-link blind steps
    permBlindCount = permBlindCount + double(permSunBlocked);  % accumulate per-link blind steps
    % Each physical ring link serves TWO terminals (terminal-2 of s1 and
    % terminal-1 of s2); throughput is counted once per terminal in use.
    capacityPerm(it) = size(permLinks,1)*2*T_P;

    % --- INTERMITTENT TERMINALS (3 & 4): three-pass bidirectional model ---

    % ---- Passes 1+2 (merged): mutual-pair identification ----------------
    % Each terminal searches BOTH adjacent planes (pPrev and pNext).
    % The hemispherical FoR naturally selects the accessible plane:
    %   A pair (p,s) <-> (q,s2) forms a canonical mutual pair
    %   T3(p,s) <-> T4(q,s2) whenever BOTH FoR conditions hold:
    %     dot(dv, +np_p) > 0   [dv in T3 hemisphere at (p,s)]
    %     dot(dv, +np_q) > 0   [-dv in T4 hemisphere at (q,s2)]
    %   where dv = pos(q,s2) - pos(p,s).
    % The complementary case (dot<0 both sides -> T4(p)<->T3(q)) is found
    % automatically when the outer loop reaches satellite q: dv flips sign,
    % both dots become positive, and the same pair is stored as T3(q)<->T4(p).
    % This canonical convention (T3-side satellite listed first) ensures
    % every physical mutual pair appears exactly once in mutualPairs.
    % FoR-mismatch pairs (dots of opposite sign) cannot form mutual links
    % and are counted in numHemiRejected.
    mutualPairs  = zeros(0, 5);  % [p, s, q, s2, dist]  T3(p,s) <-> T4(q,s2)
    numHemiRej   = 0;
    numTotalProp = 0;

    for p = 1:P
        np_p  = planeNormal(p,:)';
        pPrev = mod(p-2, P) + 1;
        pNext = mod(p,   P) + 1;

        for s = 1:Ns
            r1 = squeeze(posAll(p,s,:));

            for qAdj = [pPrev, pNext]
                np_q = planeNormal(qAdj,:)';

                for s2 = 1:Ns
                    r2   = squeeze(posAll(qAdj,s2,:));
                    dv   = r2 - r1;
                    dist = norm(dv);

                    if dist < eps || dist > maxISLrange
                        continue;                     % trivial or out of range
                    end
                    numTotalProp = numTotalProp + 1;  % in-range candidate

                    dpP = dot(dv, np_p);
                    dpQ = dot(dv, np_q);

                    if dpP > 0 && dpQ > 0
                        % Canonical form: T3(p,s) <-> T4(qAdj,s2)
                        mutualPairs(end+1,:) = [p s qAdj s2 dist]; %#ok<AGROW>
                    elseif ~(dpP < 0 && dpQ < 0)
                        % Opposite-sign dots -> FoR mismatch, no possible mutual link
                        % (dpP<0,dpQ<0 pairs are valid but found from q's iteration)
                        numHemiRej = numHemiRej + 1;
                    end
                end
            end
        end
    end

    % ---- Pass 3: greedy nearest-first assignment with Earth + Sun gates ----
    % Sort mutual pairs by distance (nearest first). Accept each pair only
    % if both terminals are still free AND the link passes Earth-occlusion
    % AND Sun-blinding checks. Once a terminal is used it cannot be
    % reassigned (at most one confirmed link per terminal).
    % Each entry in mutualPairs is canonical T3(p,s) <-> T4(q,s2), so:
    %   termFree3(p,s) and termFree4(q,s2) must both be true to confirm.
    intermLinks    = zeros(0, 4);   % [p, s, q, s2] confirmed links
    termFree3      = true(P, Ns);   % terminal-3 availability
    termFree4      = true(P, Ns);   % terminal-4 availability
    numEarthBlk    = 0;
    numSunBlk      = 0;
    numBusy        = 0;

    if ~isempty(mutualPairs)
        [~, ord]    = sort(mutualPairs(:,5), 'ascend');
        mutualPairs = mutualPairs(ord, :);
        for kk = 1:size(mutualPairs,1)
            p  = mutualPairs(kk,1);  s  = mutualPairs(kk,2);
            p2 = mutualPairs(kk,3);  s2 = mutualPairs(kk,4);
            % Gate A: terminal availability - T3 at (p,s), T4 at (p2,s2)
            if ~termFree3(p,s) || ~termFree4(p2,s2)
                numBusy = numBusy + 1;  continue;
            end
            r1 = squeeze(posAll(p, s, :));
            r2 = squeeze(posAll(p2,s2,:));
            % Gate B: Earth occlusion
            if earthBlocks(r1, r2, Re)
                numEarthBlk = numEarthBlk + 1;  continue;
            end
            % Gate C: Sun-blinding
            if sunBlinded(r1, r2, sunDir, phi, Re)
                numSunBlk = numSunBlk + 1;
                termBlindCount3(p, s)  = termBlindCount3(p, s)  + 1;  % T3 of (p,s)  blinded
                termBlindCount4(p2,s2) = termBlindCount4(p2,s2) + 1;  % T4 of (p2,s2) blinded
                continue;
            end
            % Confirm link; mark T3(p,s) and T4(p2,s2) as busy
            intermLinks(end+1,:) = [p s p2 s2]; %#ok<AGROW>
            termFree3(p, s)   = false;
            termFree4(p2,s2)  = false;
        end
    end

    numHemiRejected(it)    = numHemiRej;
    numTotalProposals(it)  = numTotalProp;
    numMutualPairs(it)     = size(mutualPairs,1);
    numEarthBlocked(it)    = numEarthBlk;
    numSunBlockedInterm(it)= numSunBlk;
    numBusyTerminal(it)    = numBusy;

    % Each confirmed link: one terminal each side -> 2*T_I aggregate capacity
    numIntermLinks(it) = size(intermLinks,1);
    capacityInterm(it) = numIntermLinks(it) * 2*T_I;
    capacityTotal(it)  = capacityPerm(it) + capacityInterm(it);

    % --- BISECTION BANDWIDTH & OVERSUBSCRIPTION RATIO -----------------
    planeCapacitySum = zeros(P,1);
    for kk = 1:size(permLinks,1)
        p = permLinks(kk,1);
        planeCapacitySum(p) = planeCapacitySum(p) + 2*T_P;
    end
    % Each confirmed inter-ring link occupies one terminal on EACH endpoint,
    % so both planes gain T_I of attached capacity.
    for kk = 1:size(intermLinks,1)
        p1 = intermLinks(kk,1);  p2 = intermLinks(kk,3);
        planeCapacitySum(p1) = planeCapacitySum(p1) + T_I;
        planeCapacitySum(p2) = planeCapacitySum(p2) + T_I;
    end

    cutCapAll   = zeros(P,1);
    idealCapAll = zeros(P,1);
    for r = 1:P
        gA = maskGroupA(r,:);
        idealCapAll(r) = sum(planeCapacitySum(gA));
        if ~isempty(intermLinks)
            % Each confirmed link crossing the cut carries 2*T_I of
            % bidirectional capacity across the boundary.
            crossing = gA(intermLinks(:,1)) ~= gA(intermLinks(:,3));
            cutCapAll(r) = 2*T_I * sum(crossing);
        end
    end

    [bbVal, rStar] = min(cutCapAll);
    bisectionBW(it) = bbVal;
    if bbVal > 0
        oversubRatio(it) = idealCapAll(rStar) / bbVal;
    else
        oversubRatio(it) = NaN;
    end

    % --- LATENCY & RELIABILITY ----------------------------------------
    % Build a weighted undirected graph of ALL currently active ISLs.
    % Edge weight = propagation delay [ms] = link_distance_km / c_light.
    % Node index for satellite (p,s): (p-1)*Ns + s  (1 ... P*Ns).
    nEdges  = size(permLinks,1) + size(intermLinks,1);
    eI = zeros(nEdges,1);  eJ = zeros(nEdges,1);  eW = zeros(nEdges,1);
    eidx = 0;

    for kk = 1:size(permLinks,1)
        p  = permLinks(kk,1); s1 = permLinks(kk,2); s2 = permLinks(kk,3);
        r1 = squeeze(posAll(p,s1,:));  r2 = squeeze(posAll(p,s2,:));
        eidx = eidx+1;
        eI(eidx) = (p-1)*Ns+s1;  eJ(eidx) = (p-1)*Ns+s2;
        eW(eidx) = norm(r1-r2) / c_light * 1e3;  % [ms]
    end
    for kk = 1:size(intermLinks,1)
        p1=intermLinks(kk,1); s1=intermLinks(kk,2);
        p2=intermLinks(kk,3); s2=intermLinks(kk,4);
        r1=squeeze(posAll(p1,s1,:)); r2=squeeze(posAll(p2,s2,:));
        eidx = eidx+1;
        eI(eidx) = (p1-1)*Ns+s1;  eJ(eidx) = (p2-1)*Ns+s2;
        eW(eidx) = norm(r1-r2) / c_light * 1e3;  % [ms]
    end

    % Fraction of maximum possible links that are active.
    % Max perm = P*Ns (ring); max inter-ring = P*Ns (one T3 per satellite).
    fracLinksActive(it) = eidx / (2*P*Ns);

    if eidx > 0
        % All-pairs shortest paths via MATLAB graph/distances (Dijkstra).
        G  = graph(eI(1:eidx), eJ(1:eidx), eW(1:eidx), N_nodes);
        D  = distances(G);          % N_nodes x N_nodes, Inf if unreachable

        % Single-hop mean delay
        meanHopLatency(it) = mean(eW(1:eidx));

        % Off-diagonal entries only
        offDiagMask = ~eye(N_nodes,'logical');
        Dvec  = D(offDiagMask);     % N_nodes*(N_nodes-1) values
        reach = isfinite(Dvec);

        fracReachable(it)  = mean(reach);
        netConnected(it)   = all(reach);
        if any(reach)
            meanE2ELatency(it) = mean(Dvec(reach));
            maxE2ELatency(it)  = max(Dvec(reach));
        else
            meanE2ELatency(it) = NaN;
            maxE2ELatency(it)  = NaN;
        end

        % Hop counts: rebuild with unit weights so distances() returns hops.
        G_hops = graph(eI(1:eidx), eJ(1:eidx), ones(eidx,1), N_nodes);
        H      = distances(G_hops);          % N_nodes x N_nodes, hop counts
        Hvec   = H(offDiagMask);
        hreach = isfinite(Hvec);
        if any(hreach)
            meanHops(it) = mean(Hvec(hreach));
            maxHops(it)  = max(Hvec(hreach));
        else
            meanHops(it) = NaN;
            maxHops(it)  = NaN;
        end
    else
        meanHopLatency(it) = NaN;
        meanE2ELatency(it) = NaN;
        maxE2ELatency(it)  = NaN;
        fracReachable(it)  = 0;
        netConnected(it)   = false;
        meanHops(it)       = NaN;
        maxHops(it)        = NaN;
    end

    % --- PLOTTING (refresh every frame) ---
    allPos = reshape(posAll, P*Ns, 3);
    litFlat = reshape(litMask, P*Ns, 1);
    set(hSat, 'XData',allPos(litFlat,1), 'YData',allPos(litFlat,2), 'ZData',allPos(litFlat,3));
    set(hSatEclipsed, 'XData',allPos(~litFlat,1), 'YData',allPos(~litFlat,2), 'ZData',allPos(~litFlat,3));

    delete(hPerm); delete(hInt); delete(hPermBlocked);

    hPerm = gobjects(size(permLinks,1),1);
    for kk = 1:size(permLinks,1)
        p = permLinks(kk,1); s1 = permLinks(kk,2); s2 = permLinks(kk,3);
        r1 = squeeze(posAll(p,s1,:)); r2 = squeeze(posAll(p,s2,:));
        hPerm(kk) = plot3(hAx,[r1(1) r2(1)],[r1(2) r2(2)],[r1(3) r2(3)], ...
                           'g-','LineWidth',1.0);
    end

    blockedPermLinks = permLinksAll(permSunBlocked,:);
    hPermBlocked = gobjects(size(blockedPermLinks,1),1);
    for kk = 1:size(blockedPermLinks,1)
        p = blockedPermLinks(kk,1); s1 = blockedPermLinks(kk,2); s2 = blockedPermLinks(kk,3);
        r1 = squeeze(posAll(p,s1,:)); r2 = squeeze(posAll(p,s2,:));
        hPermBlocked(kk) = plot3(hAx,[r1(1) r2(1)],[r1(2) r2(2)],[r1(3) r2(3)], ...
                           'Color',[1 0.6 0],'LineStyle',':','LineWidth',1.2);
    end

    hInt = gobjects(size(intermLinks,1),1);
    for kk = 1:size(intermLinks,1)
        p1=intermLinks(kk,1); s1=intermLinks(kk,2);
        p2=intermLinks(kk,3); s2=intermLinks(kk,4);
        r1 = squeeze(posAll(p1,s1,:)); r2 = squeeze(posAll(p2,s2,:));
        hInt(kk) = plot3(hAx,[r1(1) r2(1)],[r1(2) r2(2)],[r1(3) r2(3)], ...
                          'r-','LineWidth',1.2);
    end

    title(hAx, sprintf(['LEO Constellation (4-terminal)  t = %.0f s (%.1f min)  |  ' ...
                  'Eclipsed = %d/%d  |  Connected: %s\n' ...
                  'Perm links = %d (blinded = %d)  |  Inter-ring confirmed = %d  |  ' ...
                  'Capacity = %.0f Gbps  |  Bisect BW = %.0f Gbps\n' ...
                  'Hop delay = %.2f ms  |  Mean E2E = %.2f ms  |  ' ...
                  'Max E2E = %.2f ms  |  Reachable pairs = %.1f%%'], ...
                  t, t/60, numEclipsed(it), P*Ns, ...
                  mat2str(netConnected(it)), ...
                  size(permLinks,1), numSunBlockedPerm(it), numIntermLinks(it), ...
                  capacityTotal(it), bisectionBW(it), ...
                  meanHopLatency(it), meanE2ELatency(it), ...
                  maxE2ELatency(it), fracReachable(it)*100));
    drawnow limitrate;
end

%% ---------------- POST-PROCESSING PLOTS -----------
figure('Color','w');
plot(tVec/60, capacityPerm,   'g-','LineWidth',1.5); hold on;
plot(tVec/60, capacityInterm, 'r-','LineWidth',1.5);
plot(tVec/60, capacityTotal,  'k-','LineWidth',2);
xlabel('Time [min]'); ylabel('Aggregate network capacity [Gbps]');
legend('Permanent links','Intermittent terminal-links','Total','Location','best');
title('Network Capacity vs Time (4-terminal model)'); grid on;

figure('Color','w');
subplot(2,1,1);
plot(tVec/60, numTotalProposals, 'b-',  'LineWidth',1.5); hold on;
plot(tVec/60, numMutualPairs,    'm-',  'LineWidth',1.5);
plot(tVec/60, numIntermLinks,    'g-',  'LineWidth',2.0);
xlabel('Time [min]'); ylabel('Count');
legend('Total proposals (range+FoR)','Mutual pairs','Confirmed links','Location','best');
title('Inter-ring Link Pipeline vs Time'); grid on;
subplot(2,1,2);
plot(tVec/60, numHemiRejected,    'b--', 'LineWidth',1.5); hold on;
plot(tVec/60, numEarthBlocked,    'Color',[0.3 0.3 0.3],'LineWidth',1.5);
plot(tVec/60, numSunBlockedInterm,'Color',[0.6 0.2 0.8],'LineWidth',1.5);
plot(tVec/60, numBusyTerminal,    'Color',[0.1 0.6 0.6],'LineWidth',1.5);
xlabel('Time [min]'); ylabel('Count');
legend('FoR-rejected (Pass 1)','Earth-blocked (Pass 3)', ...
       'Sun-blinded (Pass 3)','Busy terminal (Pass 3)','Location','best');
title('Rejection Breakdown vs Time'); grid on;

figure('Color','w');
plot(tVec/60, numSunBlockedPerm,   'Color',[1 0.6 0],'LineWidth',1.5); hold on;
plot(tVec/60, numSunBlockedInterm, 'Color',[0.6 0.2 0.8],'LineWidth',1.5);
plot(tVec/60, numEclipsed,         'Color',[0.4 0.4 0.4],'LineWidth',1.5);
xlabel('Time [min]'); ylabel('Count');
legend('Sun-blinded permanent links','Sun-blinded intermittent cand.', ...
       'Eclipsed satellites','Location','best');
title(sprintf('Sun-Outage / Eclipse vs Time (\\phi = %.0f^\\circ, Sun lon. = %.0f^\\circ, lat. = %.0f^\\circ)', ...
      phi_deg, sunLon_deg, sunLat_deg));
grid on;

figure('Color','w');
subplot(2,1,1);
plot(tVec/60, bisectionBW, 'b-','LineWidth',1.5);
xlabel('Time [min]'); ylabel('Bisection bandwidth [Gbps]');
title('Worst-case (minimum-cut) Bisection Bandwidth vs Time'); grid on;
subplot(2,1,2);
plot(tVec/60, oversubRatio, 'Color',[0.85 0.33 0.1],'LineWidth',1.5);
xlabel('Time [min]'); ylabel('Oversubscription ratio [-]');
title('Oversubscription Ratio vs Time'); grid on;

avgBisectionBW  = mean(bisectionBW);
avgOversubRatio = mean(oversubRatio(~isnan(oversubRatio)));

% --- LATENCY PLOTS ---
figure('Color','w');
subplot(3,1,1);
plot(tVec/60, meanHopLatency, 'b-','LineWidth',1.5);
xlabel('Time [min]'); ylabel('Delay [ms]');
title('Mean Single-Hop ISL Propagation Delay'); grid on;
subplot(3,1,2);
plot(tVec/60, meanE2ELatency, 'm-','LineWidth',1.5);
xlabel('Time [min]'); ylabel('Delay [ms]');
title('Mean End-to-End Shortest-Path Propagation Delay (reachable pairs)'); grid on;
subplot(3,1,3);
plot(tVec/60, maxE2ELatency, 'r-','LineWidth',1.5);
xlabel('Time [min]'); ylabel('Delay [ms]');
title('Network Diameter (Max End-to-End Propagation Delay)'); grid on;

% --- RELIABILITY PLOTS ---
figure('Color','w');
subplot(3,1,1);
plot(tVec/60, double(netConnected), 'k-','LineWidth',1.5);
ylim([-0.05 1.05]); yticks([0 1]); yticklabels({'Partitioned','Connected'});
xlabel('Time [min]'); title('Network Connectivity vs Time'); grid on;
subplot(3,1,2);
plot(tVec/60, fracReachable*100, 'b-','LineWidth',1.5);
xlabel('Time [min]'); ylabel('Reachable pairs [%]');
title('Fraction of Reachable Satellite Pairs vs Time'); grid on;
subplot(3,1,3);
plot(tVec/60, fracLinksActive*100, 'Color',[0 0.6 0.3],'LineWidth',1.5);
xlabel('Time [min]'); ylabel('Active links [%]');
title(sprintf('Link Activity vs Time  (max = %d perm + %d inter-ring)', P*Ns, P*Ns));
grid on;

figure('Color','w');
plot(tVec/60, meanHops, 'b-','LineWidth',1.5); hold on;
plot(tVec/60, maxHops,  'r-','LineWidth',1.5);
xlabel('Time [min]'); ylabel('Hops');
legend('Mean hops (reachable pairs)','Max hops / diameter (reachable pairs)','Location','best');
title('Hop Count vs Time'); grid on;

fprintf('\n--- Bisection bandwidth / oversubscription summary ---\n');
fprintf('Average bisection bandwidth   : %.2f Gbps\n', avgBisectionBW);
fprintf('Average oversubscription ratio: %.2f\n', avgOversubRatio);
fprintf('\n--- Latency summary ---\n');
fprintf('Mean single-hop ISL delay     : %.3f ms\n', mean(meanHopLatency,'omitnan'));
fprintf('Mean end-to-end delay         : %.3f ms\n', mean(meanE2ELatency,'omitnan'));
fprintf('Mean network diameter         : %.3f ms\n', mean(maxE2ELatency,'omitnan'));
fprintf('Min  network diameter         : %.3f ms\n', min(maxE2ELatency));
fprintf('Max  network diameter         : %.3f ms\n', max(maxE2ELatency));
fprintf('\n--- Reliability summary ---\n');
fprintf('Network availability          : %.1f %%  (fraction of time fully connected)\n', ...
        mean(netConnected)*100);
fprintf('Mean reachable pairs          : %.1f %%\n', mean(fracReachable)*100);
fprintf('Mean link activity            : %.1f %%  (of max %d links)\n', ...
        mean(fracLinksActive)*100, 2*P*Ns);
fprintf('\n--- Inter-ring link pipeline summary ---\n');
fprintf('Max possible confirmed links        : %d  (= P*Ns)\n', P*Ns);
fprintf('Avg proposals made (range+FoR)      : %.1f / %d\n', ...
        mean(numTotalProposals), 2*P*Ns*Ns);
fprintf('Avg mutual pairs found              : %.1f\n', mean(numMutualPairs));
fprintf('Avg confirmed bidirectional links   : %.2f\n', mean(numIntermLinks));
fprintf('\n--- Rejection breakdown (per timestep averages) ---\n');
fprintf('Avg FoR-rejected (Pass 1)           : %.1f\n', mean(numHemiRejected));
fprintf('Avg Earth-blocked (Pass 3)          : %.1f\n', mean(numEarthBlocked));
fprintf('Avg Sun-blinded (Pass 3)            : %.1f\n', mean(numSunBlockedInterm));
fprintf('Avg busy-terminal skips (Pass 3)    : %.1f\n', mean(numBusyTerminal));
fprintf('\n--- Sun / eclipse summary (lon = %.0f deg, lat = %.0f deg, phi = %.0f deg) ---\n', ...
        sunLon_deg, sunLat_deg, phi_deg);
fprintf('Avg permanent links Sun-blinded     : %.2f / %d\n', ...
        mean(numSunBlockedPerm), P*Ns);
fprintf('Avg satellites eclipsed             : %.2f / %d\n', mean(numEclipsed), P*Ns);

%% --- Hop-count summary ---
fprintf('\n--- Hop-count summary (over reachable satellite pairs) ---\n');
fprintf('Simulation length: %d steps\n', Nt);
fprintf('Mean of mean hops per step : %.2f\n', mean(meanHops,'omitnan'));
fprintf('Mean of max  hops per step : %.2f\n', mean(maxHops, 'omitnan'));
fprintf('Overall max hops (diameter): %d\n',   max(maxHops));

%% --- Solar blinding step counts per terminal ---
% Intra-ring: each ring link kk = [p, s1, s2] blinds terminal-2 of s1
% and terminal-1 of s2 simultaneously, so permBlindCount(kk) is the
% blinding step count shared by that terminal pair.
% Inter-ring: termBlindCount3/4 count the steps in which a terminal's
% nearest mutual-pair candidate is rejected by Sun-blinding in Pass 3.
fprintf('\n--- Solar blinding step counts per terminal ---\n');
fprintf('Simulation length: %d steps\n\n', Nt);
fprintf('Intra-ring ISL terminals (terminal-1 / terminal-2, per ring link):\n');
fprintf('  Average blind steps : %.1f / %d  (%.1f %%)\n', ...
        mean(permBlindCount), Nt, mean(permBlindCount)/Nt*100);
fprintf('  Maximum blind steps : %d / %d  (%.1f %%)\n', ...
        max(permBlindCount),  Nt, max(permBlindCount)/Nt*100);
fprintf('\nInter-ring ISL terminal-3 (one per satellite, P*Ns = %d total):\n', P*Ns);
fprintf('  Average blind steps : %.1f / %d  (%.1f %%)\n', ...
        mean(termBlindCount3(:)), Nt, mean(termBlindCount3(:))/Nt*100);
fprintf('  Maximum blind steps : %d / %d  (%.1f %%)\n', ...
        max(termBlindCount3(:)),  Nt, max(termBlindCount3(:))/Nt*100);
fprintf('\nInter-ring ISL terminal-4 (one per satellite, P*Ns = %d total):\n', P*Ns);
fprintf('  Average blind steps : %.1f / %d  (%.1f %%)\n', ...
        mean(termBlindCount4(:)), Nt, mean(termBlindCount4(:))/Nt*100);
fprintf('  Maximum blind steps : %d / %d  (%.1f %%)\n', ...
        max(termBlindCount4(:)),  Nt, max(termBlindCount4(:))/Nt*100);

%% --- Solar blinding step counts per terminal ---
% Intra-ring: each ring link kk (= [p, s1, s2]) blinds terminal-2 of s1
% and terminal-1 of s2 simultaneously, so permBlindCount(kk) is the
% blinding step count shared by that terminal pair.
fprintf('\n--- Solar blinding step counts ---\n');
fprintf('Simulation length: %d steps\n\n', Nt);
fprintf('Intra-ring ISL terminals (terminal-1 / terminal-2, per ring link):\n');
fprintf('  Average blind steps : %.1f / %d  (%.1f %%)\n', ...
        mean(permBlindCount), Nt, mean(permBlindCount)/Nt*100);
fprintf('  Maximum blind steps : %d / %d  (%.1f %%)\n', ...
        max(permBlindCount),  Nt, max(permBlindCount)/Nt*100);
fprintf('\nInter-ring ISL terminal-3 (one per satellite, P*Ns = %d total):\n', P*Ns);
fprintf('  Average blind steps : %.1f / %d  (%.1f %%)\n', ...
        mean(termBlindCount3(:)), Nt, mean(termBlindCount3(:))/Nt*100);
fprintf('  Maximum blind steps : %d / %d  (%.1f %%)\n', ...
        max(termBlindCount3(:)),  Nt, max(termBlindCount3(:))/Nt*100);
fprintf('\nInter-ring ISL terminal-4 (one per satellite, P*Ns = %d total):\n', P*Ns);
fprintf('  Average blind steps : %.1f / %d  (%.1f %%)\n', ...
        mean(termBlindCount4(:)), Nt, mean(termBlindCount4(:))/Nt*100);
fprintf('  Maximum blind steps : %d / %d  (%.1f %%)\n', ...
        max(termBlindCount4(:)),  Nt, max(termBlindCount4(:))/Nt*100);

%% ---------------- HELPER FUNCTIONS -----------------
function blocked = earthBlocks(r1, r2, Re)
% Simple line-of-sight occultation check: returns true if the straight
% segment between r1 and r2 passes through the Earth sphere of radius Re.
    d  = r2 - r1;
    dd = dot(d,d);
    if dd < eps
        blocked = false; return;
    end
    tmin = -dot(r1,d)/dd;
    tmin = max(0,min(1,tmin));
    closest = r1 + tmin*d;
    blocked = norm(closest) < Re;
end

function blinded = sunBlinded(r1, r2, sunDir, phi, Re)
% An ISL is considered unavailable due to the Sun only if BOTH of the
% following hold:
%   (1) its line-of-sight forms an angle smaller than phi [rad] with the
%       Earth-Sun line (the optical/RF terminal would have to point too
%       close to the Sun), and
%   (2) BOTH endpoints are currently in sunlight (i.e. not eclipsed by
%       the Earth's shadow). A satellite in Earth's shadow cannot be
%       "blinded by the Sun" since it has no direct view of it.
    d = r2 - r1;
    nd = norm(d);
    if nd < eps
        blinded = false; return;
    end
    uLink = d / nd;
    cosAngle = abs(dot(uLink, sunDir));   % sunDir assumed unit vector
    cosAngle = min(1, max(-1, cosAngle)); % numerical safety
    angle = acos(cosAngle);               % in [0, pi/2]

    anglePass = angle < phi;
    if ~anglePass
        blinded = false; return;          % short-circuit: angle test already fails
    end

    bothSunlit = inSunlight(r1, sunDir, Re) && inSunlight(r2, sunDir, Re);
    blinded = anglePass && bothSunlit;
end

function lit = inSunlight(r, sunDir, Re)
% Simple cylindrical Earth-shadow (eclipse) model: a satellite is in
% Earth's shadow (umbra) if it is on the night side of the Earth-Sun
% line (component along -sunDir) AND its perpendicular distance from
% that line is smaller than the Earth's radius Re. Otherwise it is lit.
    along = dot(r, sunDir);               % projection onto Sun direction
    perp  = r - along*sunDir;             % component perpendicular to it
    inShadowCylinder = norm(perp) < Re;
    onNightSide = along < 0;
    lit = ~(onNightSide && inShadowCylinder);
end