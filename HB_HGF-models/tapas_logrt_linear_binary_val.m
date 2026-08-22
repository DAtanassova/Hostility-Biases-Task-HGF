function [logp, yhat, res] = tapas_logrt_linear_binary_val(r, infStates, ptrans)
% Calculates the log-probability of log-reaction times y (in units of log-ms) according to the
% linear log-RT model developed with Louise Marshall and Sven Bestmann
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2014-2016 Christoph Mathys, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Transform parameters to their native space
be0  = ptrans(1);   % intercept
be1  = ptrans(2);   % main effect: condition (hostile)
be2  = ptrans(3);   % main effect: surprise
be3  = ptrans(4);   % interaction: condition × surprise
be4  = ptrans(5);   % bernoulli variance
be5  = ptrans(6);   % inferential variance
be6  = ptrans(7);   % phasic volatility
ze   = exp(ptrans(8)); % noise

% Number of trials
n = size(infStates,1);

% Inputs
u = r.u(:,1);
% Valence
val = r.u(:,2);


% Weed irregular trials out from responses and inputs
y = r.y(:,2);
y = log(y);
y(r.irr) = [];
u(r.irr) = [];
val(r.irr) = [];


% Extract trajectories of interest from infStates
mu1hat = infStates(:,1,1);
sa1hat = infStates(:,1,2);
mu2    = infStates(:,2,3);
sa2    = infStates(:,2,4);
mu3    = infStates(:,3,3);

% Surprise
% ~~~~~~~~
m1hreg = mu1hat;
m1hreg(r.irr) = [];
poo = m1hreg.^u.*(1-m1hreg).^(1-u); % probability of observed outcome
surp = -log2(poo);

% Bernoulli variance (aka irreducible uncertainty, risk) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bernv = sa1hat;
bernv(r.irr) = [];

% Inferential variance (aka informational or estimation uncertainty, ambiguity)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level
inferv(r.irr) = [];

% Phasic volatility (aka environmental or unexpected uncertainty)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pv = tapas_sgm(mu2, 1).*(1-tapas_sgm(mu2, 1)).*exp(mu3); % transform down to 1st level
pv(r.irr) = [];

% Calculate predicted log-reaction time
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% logrt = be0 + be1.*bernv +be2.*inferv +be3.*pv;
logrt = be0 + be1.*val + be2.*surp + be3.*(val.*surp) + be4.*bernv +  be5.*inferv + be6.*pv;

% Calculate log-probabilities for non-irregular trials
% Note: 8*atan(1) == 2*pi (this is used to guard against
% errors resulting from having used pi as a variable).
reg = ~ismember(1:n,r.irr);
logp(reg) = -1/2.*log(8*atan(1).*ze) -(y-logrt).^2./(2.*ze);
yhat(reg) = logrt;
res(reg) = y-logrt;

return;
