function [rel_beam, cons] = RelatBeam(BeamVoltage, DeltaV)
% This function determines the properties of a relativistic electron beam. 
%
% Usage:
% [beam, cons] = RelatBeam(BeamVoltage, DeltaV)
% BeamVoltage: the accelerating voltage of the beam
% DeltaV: is the FWHM of the energy spread of the beam 
%         (can be zero for consider monochromatic case)
%
% beam: contains properties of the beam. See in the function.
% cons: contains the physical constants used in the calculations 
% (all in SI units).
%
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2016
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************


rel_beam.voltage = BeamVoltage;
rel_beam.del_V = DeltaV;

% Constants:
cons.m0    = 9.109382e-31;
cons.c     = 299792458;
cons.ec    = 1.602176e-19;
cons.h     = 6.626069e-34;
cons.mu_0  = 4*pi*1e-7;
cons.eps_0 = 1/(cons.c^2*cons.mu_0);

% Physical constant for relativisitic corrections:
cons.epsilon = cons.ec/(2*cons.m0*cons.c^2);
epsilon = cons.epsilon;

% epsilon = 0.97845e-6; % should be this!
% and where...
% Vr = V*(1+epsilon*V) and 
% energy spread Vr = DeltaV * ( 1 + 2*epsilon*V)
% or in fact more accurately, by just doing an expansion of DeltaVrel below:
% energy spread Vr = DeltaV * ( 1 + 2*epsilon*V + epsilon*DeltaV),
% But of course generally DeltaV very small compared to V if we are worried
% about relativistic effects.
% Also,
% p = hbar/k = h/lambda
% eV = p^2/(2m) = h^2/(2m*lambda^2);
% lambda = h/sqrt(2eVm)
cons.vol_lambda = cons.h/sqrt(2*cons.m0*cons.ec);
% should be ~1.22639e-9

% The relativistic beam voltage, at the centre energy.
BeamVoltageRel_cen = BeamVoltage + epsilon*BeamVoltage^2;

% The relativistic beam voltage, at the offset energy.
BeamVoltage_off = BeamVoltage + DeltaV;
BeamVoltageRel_off = BeamVoltage_off + epsilon*BeamVoltage_off.^2;
% and the relativisitic delta, this should be the same as above?
rel_beam.BeamVoltageRel_cen = BeamVoltageRel_cen;
rel_beam.BeamVoltageRel_off = BeamVoltageRel_off;
rel_beam.DeltaVRel = BeamVoltageRel_off - BeamVoltageRel_cen;

% We actually should be working with the lambda for the actual energy in
% question, not the centre voltage - not sure how many people actually do
% this though...
rel_beam.BeamLambdaRel = cons.vol_lambda./sqrt(BeamVoltageRel_off);
rel_beam.KRel = 2*pi./rel_beam.BeamLambdaRel;

rel_beam.mRel = (cons.ec*BeamVoltage/cons.c^2) + cons.m0;

rel_beam.beta = (1-(((cons.ec*BeamVoltage/(cons.m0*cons.c^2))+1)^-2))^0.5;

rel_beam.vRel = rel_beam.beta * cons.c;

rel_beam.sigma  = ...
    (2 * pi * rel_beam.mRel * cons.ec* rel_beam.BeamLambdaRel / (cons.h^2));

rel_beam.p_rel = rel_beam.vRel.*rel_beam.mRel;
% The velocity that would give a rest mass electron the correct momentum.
% can be handy in modelling magnetic lenses without too much fuss
rel_beam.equiv_vel = rel_beam.p_rel/cons.m0;

