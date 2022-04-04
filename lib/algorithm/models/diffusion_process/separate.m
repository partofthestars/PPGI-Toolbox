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
%	separate.m:
%
% Description:
%
%   Implements the Kalman based separation of frequencies based upon the diffusion process algorithm.
%
% References:
%
%	Christian S. Pilz, Jarek Krajewski, Vladimir Blazek.
%   On the Diffusion Process for Heart Rate Estimation from Face Videos under Realistic Conditions.
%   Pattern Recognition: 39th German Conference, GCPR 2017, Basel, Switzerland.
%   Proceedings (Lecture Notes in Computer Science), Springer, 2017
%

function [S,CS,QCS,SMM,SPP,FMM,FPP] = separate(Y,dt,FF,nharm,BF,BQ,BL,BH,R,qr)
% - Separate periodic signal from other signals
%
% Syntax:
%   [S,CS,QCS,SMM,SPP,FMM,FPP] = separate(Y,dt,FF,nharm,BF,BQ,BL,BH,R,qr)
%
% In:
%   Y  - Signal as 1xD vector
%   dt - Sampling period in seconds (e.g. 0.1)
%   FF - Fundamental frequencies as 1xD vector in BPM (from the IMM)
%   nharm - Number of harmonics to be estimated (e.g. 3)
%   BF - Feedback matrix for the bias model (e.g. [0 1; 0 0])
%   BQ - Process spectral density for the bias model (e.g. 0.01)
%   BL - Noise multiplier matrix for the bias model (e.g. [0;1])
%   BH - Measurement matrix for the bias model (e.g. [1 0])
%   R  - Measurement variance (e.g. 0.1^2)
%   qr - Resonator's process noise spectral density (e.g. 0.1)
%
% Out:
%   S   - Cleaned signal as 1xD vector
%   CS  - Fundamental signal and its harmonics as (nharm)xD vector
%   QCS - Quadrature periodic signals as (nharm)xD vector
%   SMM - The smoother means NxD
%   SPP - The smoother covariances NxNxD
%   FMM - The filter means NxD
%   FPP - The filter covariances NxNxD
%
% Description:
%   Separate given signal into true signal and periodic
%   signal by using the known frequency of the periodic signal.
%   See TRACK for the model used.

    % Number of harmonics (including fundamental)
    N = nharm;
    
    % Prior means and covariances
    m0 = zeros(2*N+size(BF,1),1);  % XXX: Could be function parameter
    P0 = 10*eye(size(m0,1));       % XXX: Could be function parameter

    % Mean to measurement at step 1
    m0(1) = Y(1);
    
    %
    % Form the constant matrices
    %
    F  = zeros(size(m0,1));
    H  = zeros(1,size(m0,1));
    L  = zeros(size(m0,1),N+size(BQ,1));
    Qc = zeros(N+size(BQ,1));
    
    for j=1:N
        i1 = 1+2*(j-1);
        i2 = 2+2*(j-1);
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
  
    
%% Run filter
    
    m = m0;
    P = P0;

    FMM = zeros(size(m,1),length(Y));
    FPP = zeros(size(P,1),size(P,2),length(Y));
    
    for k=1:length(Y)
        
        %
        % Update the resonators
        %
        w  = 2*pi*FF(k)/60;
        for j=1:N
            i1 = 1+2*(j-1);
            i2 = 2+2*(j-1);
            F(i1:i2,i1:i2) = [0 j*w; -j*w 0];
        end
        
        %
        % Discretize
        %
        [A,Q] = lti_disc(F,L,Qc,dt);

        %
        % Estimate with KF
        %
        [m,P] = kf_predict(m,P,A,Q);
        [m,P] = kf_update(m,P,Y(k),H,R);
    
        FMM(:,k) = m;
        FPP(:,:,k) = P;
    end
    
    
%% Run smoother

    % Allocate and set initial
    SMM = zeros(size(FMM));
    SPP = zeros(size(FPP));
    SMM(:,end)   = m;
    SPP(:,:,end) = P;
            
    for k=length(Y)-1:-1:1
        
        %
        % Update each resonator
        %
        w  = 2*pi*FF(k+1)/60;
        for j=1:N
            i1 = 1+2*(j-1);
            i2 = 2+2*(j-1);
            F(i1:i2,i1:i2) = [0 j*w; -j*w 0];
        end
        
        %
        % Discretize
        %
        [A,Q] = lti_disc(F,L,Qc,dt);
        
        % Smoother
        m_pred = A * FMM(:,k);
        P_pred = A * FPP(:,:,k) * A' + Q;
        D  = FPP(:,:,k) * A' / P_pred;
        m  = FMM(:,k) + D * (m - m_pred);
        P  = FPP(:,:,k) + D * (P - P_pred) * D';
        
        % Store results
        SMM(:,k)   = m;
        SPP(:,:,k) = P;
        
    end
    
    %
    % Extract the signals
    %
    S   = SMM(end-size(BF,1)+1,:);
    CS  = SMM(1:2:(2*nharm),:);
    QCS = SMM(2:2:(2*nharm+1),:);
 
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

 function [x,P] = kf_predict(x,P,A,Q,B,u)
% KF_PREDICT - Perform Kalman Filter prediction step
%
% Syntax:
%   [X,P] = KF_PREDICT(X,P,A,Q,B,U)
%
% In:
%   X - Nx1 mean state estimate of previous step
%   P - NxN state covariance of previous step
%   A - Transition matrix of discrete model (optional, default identity)
%   Q - Process noise of discrete model     (optional, default zero)
%   B - Input effect matrix                 (optional, default identity)
%   U - Constant input                      (optional, default empty)
%
% Out:
%   X - Predicted state mean
%   P - Predicted state covariance
%   
% Description:
%   Perform Kalman Filter prediction step. The model is
%
%     x[k] = A*x[k-1] + B*u[k-1] + q,  q ~ N(0,Q).
% 
%   The predicted state is distributed as follows:
%   
%     p(x[k] | x[k-1]) = N(x[k] | A*x[k-1] + B*u[k-1], Q[k-1])
%
%   The predicted mean x-[k] and covariance P-[k] are calculated
%   with the following equations:
%
%     m-[k] = A*x[k-1] + B*u[k-1]
%     P-[k] = A*P[k-1]*A' + Q.
%
%   If there is no input u present then the first equation reduces to
%     m-[k] = A*x[k-1]
%

  %
  % Check arguments
  %
  if nargin < 3
    A = [];
  end
  if nargin < 4
    Q = [];
  end
  if nargin < 5
    B = [];
  end
  if nargin < 6
    u = [];
  end
  
  %
  % Apply defaults
  %
  if isempty(A)
    A = eye(size(x,1));
  end
  if isempty(Q)
    Q = zeros(size(x,1));
  end
  if isempty(B) & ~isempty(u)
    B = eye(size(x,1),size(u,1));
  end

  %
  % Perform prediction
  %
  if isempty(u)
    x = A * x;
    P = A * P * A' + Q;
  else
    x = A * x + B * u;
    P = A * P * A' + Q;
  end

 function [X,P,K,IM,IS,LH] = kf_update(X,P,y,H,R)
% KF_UPDATE - Kalman Filter update step
%
% Syntax:
%   [X,P,K,IM,IS,LH] = KF_UPDATE(X,P,Y,H,R)
%
% In:
%   X - Nx1 mean state estimate after prediction step
%   P - NxN state covariance after prediction step
%   Y - Dx1 measurement vector.
%   H - Measurement matrix.
%   R - Measurement noise covariance.
%
% Out:
%   X  - Updated state mean
%   P  - Updated state covariance
%   K  - Computed Kalman gain
%   IM - Mean of predictive distribution of Y
%   IS - Covariance or predictive mean of Y
%   LH - Predictive probability (likelihood) of measurement.
%   
% Description:
%   Kalman filter measurement update step. Kalman Filter
%   model is
%
%     x[k] = A*x[k-1] + B*u[k-1] + q,  q ~ N(0,Q)
%     y[k] = H*x[k]   + r,             r ~ N(0,R)
%
%   Prediction step of Kalman filter computes predicted
%   mean m-[k] and covariance P-[k] of state:
%
%     p(x[k] | y[1:k-1]) = N(x[k] | m-[k], P-[k])
%
%   See for instance KF_PREDICT how m-[k] and P-[k] are
%   calculated. 
%
%   Update step computes the posterior mean m[k] and
%   covariance P[k]  of state given new measurement:
%
%     p(x[k] |�y[1:k]) = N(x[k] |�m[k], P[k])
%
%   Innovation distribution is defined as
%
%     p(y[k] |�y[1:k-1]) = N(y[k] | IM[k], IS[k])
%   
%   Updated mean x[k] and covarience P[k] are given by
%   the following equations (not the only possible ones):
%
%     v[k] = y[k] - H[k]*m-[k]
%     S[k] = H[k]*P-[k]*H[k]' + R[k]
%     K[k] = P-[k]*H[k]'*[S[k]]^(-1) 
%     m[k] = m-[k] + K[k]*v[k]
%     P[k] = P-[k] - K[k]*S[k]*K[k]'
%
% Example:
%   m = m0;
%   P = P0;
%   M = m0;
%   for i=1:size(Y,2)
%     [m,P] = kf_predict(m,P,A,Q);
%     [m,P] = kf_update(m,P,Y(:,i),H,R);
%     M = [M m];
%   end
%

  %
  % Check which arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end

  %
  % update step
  %
  IM = H*X;
  IS = (R + H*P*H');
  K = P*H'/IS;
  X = X + K * (y-IM);
  P = P - K*IS*K';
  if nargout > 5
    LH = gauss_pdf(y,IM,IS);
  end

  function [P,E] = gauss_pdf(X,M,S)
%GAUSS_PDF  Multivariate Gaussian PDF
%
% Syntax:
%   [P,E] = GAUSS_PDF(X,M,S)
%
% In:
%   X - Dx1 value or N values as DxN matrix
%   M - Dx1 mean of distibution or N values as DxN matrix.
%   S - DxD covariance matrix
%
% Out:
%   P - Probability of X. 
%   E - Negative logarithm of P
%   
% Description:
%   Calculate values of PDF (Probability Density
%   Function) of multivariate Gaussian distribution
%
%    N(X | M, S)
%
%   Function returns probability of X in PDF. If multiple
%   X's or M's are given (as multiple columns), function
%   returns probabilities for each of them. X's and M's are
%   repeated to match each other, S must be the same for all.
%

  if size(M,2) == 1
    DX = X-repmat(M,1,size(X,2));  
    E = 0.5*sum(DX.*(S\DX),1);
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
  elseif size(X,2) == 1
    DX = repmat(X,1,size(M,2))-M;  
    E = 0.5*sum(DX.*(S\DX),1);
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
  else
    DX = X-M;  
    E = 0.5*DX'*(S\DX);
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
  end
