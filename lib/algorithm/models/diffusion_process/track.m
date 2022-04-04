%########################################################################
%
%	- PPGI Toolbox - 
%   A MATLAB toolbox for Photoplethysmography Imaging (PPGI)
%
% Author   : Christian S. Pilz
% Company  : The Nature of Space of Time
% Date     : 07.05.2019
%
% Contact  : cpi@partofthestars.com
% Web Page : www.partofthestars.com
%
% Version  : beta0.1
%
%########################################################################
%
%	track.m:
%
% Description:
%
%   Implements the IMM based prediction of frequencies based upon the diffusion process algorithm.
%
% References:
%
%	Christian S. Pilz, Jarek Krajewski, Vladimir Blazek.
%   On the Diffusion Process for Heart Rate Estimation from Face Videos under Realistic Conditions.
%   Pattern Recognition: 39th German Conference, GCPR 2017, Basel, Switzerland.
%   Proceedings (Lecture Notes in Computer Science), Springer, 2017
%

function [FF,MM,WW] = track(Y,dt,freqlist,nharm,BF,BQ,BL,BH,R,qr,ptrans,poverall)
% - Track periodic signal with IMM
%
% Syntax:
%   [FF,MM,WW] = track(Y,dt,freqlist,nharm,BF,BQ,BL,BH,R,qr,ptrans)
%
% In:
%   Y  - Signal as Dx1 or 1xD vector
%   dt - Sampling period in seconds (e.g. 0.1)
%   freqlist - List of candidate frequencies in beats-per-min (e.g. 65:75).
%   nharm - Number of harmonics to be estimated (e.g. 3)
%   BF - Feedback matrix for the bias model (e.g. [0 1; 0 0])
%   BQ - Process spectral density for the bias model (e.g. 0.01)
%   BL - Noise multiplier matrix for the bias model (e.g. [0;1])
%   BH - Measurement matrix for the bias model (e.g. [1 0])
%   R  - Measurement variance (e.g. 0.1^2)
%   qr - Resonator's process noise spectral density (e.g. 0.01)
%   ptrans - Transition probability for model change (e.g. 0.001)
%   poverall - Overall transition jump probability (e.g. 0.00)
%
% Out:
%   FF - MMSE estimates of fundamental frequencies as 1xD vector in BPM
%   MM - Filtered means of state components as (2*nharm*dim(bias))xD matrix
%   WW - Smoothed posterior probabilities for models as NxD matrix
%
% Description:
%   Tracks cardiac signal, which is buried in signal Y using the
%   following kind of model:
%
%     Y(k) = xb(k) + xr1(k) + ... + xrn(k) + e(k)
%
%   where e(k) is the measurement noise with variance R and
%
%   - xb(k) is the bias signal, which is modeled the process
%
%        X(k) = BA X(k-1) + w_k,  w_k ~ N(0,BQ)
%       xb(k) = BH X(k),
%
%   - xrj(k) are resonators (nharm of them), each having frequency j*f0,
%     where f0 is the fundamental frequency. The model can be written as
%
%        dRj/dt = [0 2*pi*j*f0; -2*pi*j*f0 0] Rj + wj(t),
%        xrj(k) = [1 0] Rj
%
%     where each wj has spectral density qr.
%
%   - The possible frequencies f0(j) given in freqlist form a HMM with
%     transition probabilities P(j+1|j) = P(j-1|j) = ptrans.
%
%   - The estimation is based on the Interacting Multiple Models
%     (IMM) algorithm, which used Kalman filters for tracking
%     each of the modes and mixing between them.
%

% Copyright (C) 2010 Simo Särkkä and Arno Solin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

%% 

    % Report what we are doing
    fprintf('%-36s:%+35s\n','IMM Tracking of frequencies','Preparing..')

    % Change frequencies to rad/s
    ww = 2*pi*freqlist/60; % To rads/s

    %
    % Compute transition probability matrix
    %
    NM = length(ww);
    pp = 1-ptrans;
    PIJ = poverall*ones(NM) + (pp-NM*poverall)*eye(NM) + ...
                       (1-pp)/2*diag(ones(NM-1,1),1) + ...
                       (1-pp)/2*diag(ones(NM-1,1),-1);
    PIJ(1,1)   = PIJ(1,1) + (1-pp)/2;
    PIJ(NM,NM) = PIJ(NM,NM) + (1-pp)/2;
    
    % Check if the PIJ transition matrix is consistent
    if (sum(PIJ)-NM) > 1e-10 | max(abs(PIJ)) > 1, 
        warning(['The transition matrix Pij might be inconsistent. ' ...
                 'Check your transition probabilities.'])
    end
    
    %
    % Compute the IMM parameters (ML,PL,AL,QL,HL,RL)
    %
    N = nharm;

    % Prior means and covariances
    m0 = zeros(2*N+size(BF,1),1);
    P0 = 100*eye(size(m0,1));

    RL = zeros(1,1,length(ww));
    HL = zeros(1,size(m0,1),length(ww));
    AL = zeros(size(m0,1),size(m0,1),length(ww));
    QL = zeros(size(m0,1),size(m0,1),length(ww));

    ML = repmat(m0,[1 length(ww)]);
    PL = repmat(P0,[1 1 length(ww)]);

    for i=1:length(ww)
        H  = zeros(1,size(m0,1));
        F  = zeros(size(m0,1));
        L  = zeros(size(m0,1),N+size(BQ,1));
        Qc = zeros(N+size(BQ,1));
        
        %
        % Form the resonators
        %
        for j=1:N
            i1 = 1+2*(j-1);
            i2 = 2+2*(j-1);
            F(i1:i2,i1:i2) = [0 j*ww(i); -j*ww(i) 0];
            L(i2,j) = 1;
            Qc(j,j) = qr;
            H(1,i1) = 1;
        end
        
        %
        % Form the bias model
        %
        i1 = 1+2*N;
        i2 = 2+2*N;
        j1 = N+2-size(BQ,1);
        j2 = N+1;
        F(i1:i2,i1:i2)  = BF;
        L(i1:i2,j1:j2)  = BL;
        H(1,i1:i2)      = BH;    
        Qc(j1:j2,j1:j2) = BQ;
 
        %
        % Discretize
        %
        [AL(:,:,i),QL(:,:,i)] = lti_disc(F,L,Qc,dt);
        HL(:,:,i) = H;    
        RL(:,:,i) = R;
        
    end

    % Show waitbar
    %handle = waitbar(0,'Please wait...');  
    loopstart = tic;
    
    %
    % Estimate with IMM
    %
    WW = zeros(NM,length(Y));
    W  = ones(1,NM)/NM;

    FF = zeros(1,length(Y));
    MM = zeros(size(m0,1),length(Y));
    for k=1:length(Y)
        %
        % IMM prediction
        %
        [ML,PL,W] = imm_predict0(ML,PL,W,PIJ,AL,QL);
        
        %
        % IMM update
        %
        [ML,PL,W,K,IM,IS,LH] = imm_update0(ML,PL,W,Y(k),HL,RL);
        WW(:,k) = W;
        FF(k) = sum(W .* ww) / sum(W);
        MM(:,k) = imm_estimate0(ML,PL,W);
        
        % Update waitbar and show time remaining
        if rem(k,40)==0 
          secondsleft = min(toc(loopstart),Inf)/k*(length(Y)-k);
          %waitbar(k/length(Y),handle,sprintf('Estimating with IMM\nTime left: %.0f min %.0f s.', ...
          %  floor(secondsleft/60),rem(secondsleft,60)))
          fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
              sprintf('%.0f min %.0f s', ...
                floor(secondsleft/60),rem(secondsleft,60)))
        end    
    end
    
    % Dump
    % save FF_1.mat FF

    %
    % Smooth the mode probilities
    %
    thr = eps;
    SWW = WW;
    SFF = FF;
    for k=length(Y)-1:-1:1
        PW = PIJ*WW(:,k);   % p(xk+1|Y1:k)
        PW(PW < thr) = thr;
        SW = SWW(:,k+1)./PW; % p(xk+1|Y1:T) / p(xk+1|Y1:k)
        SW = PIJ'*SW;      % int p(xk+1|xk) p(xk+1|Y1:T) / p(xk+1|Y1:k) dxk+1
        SW = SW .* WW(:,k); % * p(xk|y1:k)
        SWW(:,k) = SW / sum(SW);
        SFF(k) = sum(SW' .* ww) / sum(SW);
    end

    % Dump
    % save FF_2.mat FF SWW SFF PW thr PIJ SW ww WW

    % Get rid of the waitbar
    %close(handle)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
            'Done.')
         
    % Check if we have NaN values
    if sum(isnan(SFF)) == 0
      FF = 60*SFF/2/pi;
    else
      FF = 60*FF/2/pi;
      warning('NaN in IMM smoother, using filter result instead.')
    end
	
function [A,Q] = lti_disc(F,L,Q,dt)
% LTI_DISC - Discretize LTI ODE with Gaussian Noise
%
% Syntax:
%   [A,Q] = lti_disc(F,L,Qc,dt)
%
% In:
%   F  - NxN Feedback matrix
%   L  - NxL Noise effect matrix        (optional, default identity)
%   Qc - LxL Diagonal Spectral Density  (optional, default zeros)
%   dt - Time Step                      (optional, default 1)
%
% Out:
%   A - Transition matrix
%   Q - Discrete Process Covariance
%
% Description:
%   Discretize LTI ODE with Gaussian Noise. The original
%   ODE model is in form
%
%     dx/dt = F x + L w,  w ~ N(0,Qc)
%
%   Result of discretization is the model
%
%     x[k] = A x[k-1] + q, q ~ N(0,Q)
%
%   Which can be used for integrating the model
%   exactly over time steps, which are multiples
%   of dt.
%
% History:
%   11.01.2003  Covariance propagation by matrix fractions
%   20.11.2002  The first official version.

% Copyright (C) 2002, 2003 Simo Särkkä
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

%%
  %
  % Check number of arguments
  %
  if nargin < 1
    error('Too few arguments');
  end
  if nargin < 2
    L = [];
  end
  if nargin < 3
    Q = [];
  end
  if nargin < 4
    dt = [];
  end

  if isempty(L)
    L = eye(size(F,1));
  end
  if isempty(Q)
    Q = zeros(size(F,1),size(F,1));
  end
  if isempty(dt)
    dt = 1;
  end

  %
  % Closed form integration of transition matrix
  %
  A = expm(F*dt);

  %
  % Closed form integration of covariance
  % by matrix fraction decomposition
  %
  n   = size(F,1);
  Phi = [F L*Q*L'; zeros(n,n) -F'];
  AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
  Q   = AB(1:n,:)/AB((n+1):(2*n),:);

  
function [M,P,W] = imm_predict0(M,P,W,T,A,Q,B,U)
% IMM_PREDICT0 - Perform IMM prediction step
%
% Syntax:
%   [M,P,W] = IMM_PREDICT0(M,P,W,T,A,Q,B,U)
%
% Author:
%   Simo S�rkk�, 2003
%
% In:
%   M - DxN mean state estimates of previous step
%   P - DxDxN state covariances from previous step
%   W - 1xN weights (posterior probabilities of modes)
%   T - Mode transition probabilities T(i,j) = p(i|j)
%   A - DxD transition matrix of discrete model or DxDxN
%       matrix containing separate matrices for each mode
%       (optional, default identity)
%   Q - DxD process noise covariance of discrete model or DxDxN
%       matrix containing separate covariances for each mode
%       (optional, default zero matrix)
%   B - DxS input effect matrix or DxSxN matrix containing separate
%       input effect matrices for each mode (optional, default identity)
%   U - Sx1 constant input or SxN matrix containing separate
%       inputs for each mode                (optional, default empty)
%
% Out:
%   X  - Predicted state means as DxN matrix
%   P  - Predicted state covariance as DxDxN matrix
%   W  - Predicted weights as 1xN matrix
%   
% Description:
%   Perform Interacting Multiple Model (IMM)
%   prediction step. The model is
%
%     p(i|j) = T_{ij}
%     x_i[k] = A_i*x_i[k-1] + B*u[k-1] + q_i,  q_i ~ N(0,Q_i)
%     i = 1,...,N
%
% See also:
%   IMM_UPDATE, KF_PREDICT, LTI_DISC, EKF_PREDICT, EKF_UPDATE
%
% History:
%   18.02.2003  The first official version.

% Copyright (C) 2003 Simo Särkkä
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

%%
  %
  % Check arguments
  %
  if nargin < 5
    A = [];
  end
  if nargin < 6
    Q = [];
  end
  if nargin < 7
    B = [];
  end
  if nargin < 8
    U = [];
  end
  
  %
  % Apply defaults
  %
  if isempty(A)
    A = eye(size(M,1));
  end
  if isempty(Q)
    Q = zeros(size(M,1));
  end
  if isempty(B) && ~isempty(U)
    B = eye(size(M,1),size(U,1));
  end
  if size(A,3)==1
    A = repmat(A,[1 1 size(M,2)]);
  end
  if size(Q,3)==1
    Q = repmat(Q,[1 1 size(M,2)]);
  end

  %
  % Calculate mixing probabilities
  %
  mu_ij = T.*repmat(W',1,size(W,2));
  CJ = sum(mu_ij,1);
  % CJ(find(CJ==0)) = eps;
  CJ(CJ==0) = eps;   % Equivalent to find() above
  mu_ij = mu_ij ./ repmat(CJ,size(mu_ij,1),1);
  
  %
  % Mixing
  %
  M0 = zeros(size(M));
  P0 = zeros(size(P));

  for j=1:size(M,2)
    for i=1:size(M,2)
        if mu_ij(i,j) > eps %!!!!
            M0(:,j) = M0(:,j) + M(:,i)*mu_ij(i,j);
        end
    end
    for i=1:size(M,2)
        if mu_ij(i,j) > eps %!!!!
            C = P(:,:,i) + (M(:,i) - M0(:,j))*(M(:,i) - M0(:,j))';
            P0(:,:,j) = P0(:,:,j) + mu_ij(i,j) * C;
        end
    end
  end

  %
  % Mode matched prediction
  %
  for i=1:size(M,2)
    if isempty(U)
      M(:,i) = A(:,:,i) * M0(:,i);
    else
      M(:,i) = A(:,:,i) * M0(:,i) + B(:,:,i) * U(:,i);
    end
    P(:,:,i) = A(:,:,i) * P0(:,:,i) * A(:,:,i)' + Q(:,:,i);
  end

  %
  % Weights are actually the
  % normalization constants
  %
  W = CJ ./ sum(CJ);


  function [M,P,W,K,IM,IS,LH] = imm_update0(M,P,W,y,H,R)
% IMM_UPDATE0 - IMM update step
%
% Syntax:
%   [M,P,W,K,IM,IS,LH] = IMM_UPDATE0(M,P,W,Y,H,R)
%
% In:
%   M - Dx1xN mean state estimates
%   P - DxDxN state covariances
%   W - 1xN vector of mode weights
%   Y - Ex1 measurement vector.
%   H - ExD measurement matrix or ExDxN matrix
%       of mode conditional matrices
%   R - ExE measurement noise covariance or ExExN
%       matrix of mode conditional covariances
%
% Out:
%   X  - Updated state means
%   P  - Updated state covariances
%   W  - Updated mode weights
%   K  - Computed Kalman gains
%   IM - Means of predictive distributions of Y.
%   IS - Covariances or predictive means of Y.
%   LH - Predictive probability (likelihood) of measurement.
%   
% Description:
%   Perform Interactive Multiple Model (IMM) update step.
%   The IMM measurement model is:
%
%     y[k] = H_i*x_i[k] + r_i,  r_i ~ N(0,R_i)
%     i = 1,...,N
%
%   Predictive measurement distribution is defined as
%
%     p(y[k] |�y[1:k-1]) = sum w_i N(y[k] | IM_i[k], IS_i[k])
%
% See also:
%   IMM_PREDICT, KF_UPDATE
%
% History:
%   18.02.2003  The first official version.

% Copyright (C) 2003 Simo Särkkä
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%


  %
  % Check which arguments are there
  %
  if nargin < 6
    error('Too few arguments');
  end

  %
  % Duplicate singleton measurement parameters
  %
  if size(H,3)==1
    H = repmat(H,[1 1 size(M,2)]);
  end
  if size(R,3)==1
    R = repmat(R,[1 1 size(M,2)]);
  end

  %
  % Mode conditional update steps
  %
  IM = zeros(size(y,1),size(M,2));
  IS = zeros(size(y,1),size(y,1),size(M,2));
  K  = zeros(size(M,1),size(y,1),size(M,2));
  LH = zeros(1,size(M,2));

  for i=1:size(M,2)
    IM(:,i)   = H(:,:,i) * M(:,i);
    IS(:,:,i) = R(:,:,i) + H(:,:,i) * P(:,:,i) * H(:,:,i)';
    K(:,:,i)  = P(:,:,i) * H(:,:,i)' / IS(:,:,i);
    M(:,i)    = M(:,i) + K(:,:,i) * (y - IM(:,i));
    P(:,:,i)  = P(:,:,i) - K(:,:,i) * IS(:,:,i) * K(:,:,i)';
    %LH(1,i)   = gauss_pdf(y,IM(:,i),IS(:,:,i));
    DX = y-IM(:,i);  
    E = 0.5*DX'*(IS(:,:,i)\DX);
%    E = E + 0.5 * size(M,1) * log(2*pi) + 0.5 * log(det(IS(:,:,i)));
    E = E + 0.5 * log(det(IS(:,:,i)));
    LH(1,i) = exp(-E);
  end

  %
  % Mode probability update
  %
  W  = W .* LH;
  LH = sum(W);
  W  = W ./ LH;

  function [m,C] = imm_estimate0(M,P,W)
% IMM_ESTIMATE0 - Calculate mean and covariance of IMM mixture
%
% Syntax:
%   [m,C] = IMM_ESTIMATE0(M,P,W)
%
% Author:
%   Simo S�rkk�, 2003
%
% In:
%   M - Dx1xN mean state estimates
%   P - DxDxN state covariances
%   W - 1xN vector of mode weights
%
% Out:
%   m - State estimate mean
%   C - State estimate covariance
%   
% Description:
%   Calculate mean and covariance of IMM mixture.
%   This is just a helper routine to ease the interpretation
%   of mixture Gaussian output of IMM filter.
%
% See also:
%   IMM_PREDICT, IMM_UPDATE
%
% History:
%   18.02.2003  The first official version.

% Copyright (C) 2003 Simo Särkkä
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

%%
  m = zeros(size(M,1),1);
  for i=1:size(M,2)
    m = m + W(i) * M(:,i);
  end
  C = zeros(size(M,1),size(M,1));
  for i=1:size(M,2)
    C = C + W(i) * (P(:,:,i) + (M(:,i) - m) * (M(:,i) - m)');
  end

