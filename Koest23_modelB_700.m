% Function to derive POC [mg m^-3] using particulate backscattering at 700
% nm [m^-1] and chlorophyll-a concentration [mg m^-3] based on global
% dataset described in Koestner et al. 2022.

% Koestner D, Stramski D and Reynolds RA (2022) A Multivariable Empirical
% Algorithm for Estimating Particulate Organic Carbon Concentration in
% Marine Environments From Optical Backscattering and Chlorophyll-a
% Measurements. Front. Mar. Sci. 9:941950. doi: 10.3389/fmars.2022.941950

% Inputs: column vector of particulate backscattering at 700 nm (bbp
% [m^-1]) and corresponding chlorophyl-a concentration (chla [mg m^-3]); an
% alignment factor for difference between ECO and HS sensors (0.9 is
% probably appropriate but leave as 1 to be conservative); choice of low
% POC bias correction function (linear or power); alpha for prediction intervals 

% Outputs: POC and prediction interval of POC in [mg m^-3]

% Koestner et al. 2023 version uses all data with some slightly different
% restrictions on data, with new adjusted TChla from ANT26 and additional
% deep water BIOSOPE. It also has a power function bias correction (no
% negative data)

% NOTE PAPER IS UNDER REVISION AS OF August 11 2023 BUT DO NOT EXPECT
% CHANGES TO COEFFICIENTS

% V1: 2023.01.24
function [POC,PI]=Koest23_modelB_700(bbp,chla,alignment,bias,alpha)

if nargin==2
    alignment=1;
    bias='power';
    alpha=0.125; % two-tailed 75% confidence
elseif nargin==3
    bias='power';
    alpha=0.125;
end

bbp=bbp*alignment; %alignment adjustment for use with different bbp sensor

comp=chla./bbp;
comp(comp<=1)=1; %fix any comp values less than 1 to be 1, usually only an issue if chla=0


% model B coefficients
% k_modelB=[89.4229308852505; 0.188074098663738; 0.759083352712280; 0.193397796405782]; %2022 version
% bias_modelB_coef=[1.63582494446266; -21.2269735484462]; %2022 version

% % % % model coefficients
k_modelB=[52.8187501942431; 0.135288289126603; 0.884851394851513; 0.226810214797258];

if bias=='power'
    bias_modelB_coef=[1.46918996207386, -0.734453171035830];
    bias_lim=36.8;
elseif bias=='linear'
    bias_modelB_coef=[1.42263442605374, -14.7264359822160];
    bias_lim=34.8;
end

% % % % prediction interval stuff
S=[0.0248767643354040	0.00967335656242312	-0.0105363562351115	-0.00404078515530673
0.00967335656242312	0.00405822189121109	-0.00402605221279587	-0.00166729112989270
-0.0105363562351115	-0.00402605221279587	0.00483637706215940	0.00180593541611042
-0.00404078515530673	-0.00166729112989270	0.00180593541611042	0.000727608368784931]; %covariance matrix of coefficients

mse=0.0311478519483073; % mean square error from regression model

df=407-4-2; % take away 6 for model coefficients
% t=1.96556145883247; % 95% confidence two-tail, inverse of Student's t cumulative distribution function
t=tinv(1-alpha,403);


% calculate POC from multivariable model
pocB=k_modelB(1)*(bbp.^k_modelB(2)).*(comp.^k_modelB(3)).*(comp.^(k_modelB(4)*log10(bbp)));

% continious bias correction for modeled POC < bias_lim
if bias=='power'
    pocB(pocB<bias_lim)=(pocB(pocB<bias_lim).^bias_modelB_coef(1))*10^bias_modelB_coef(2);
elseif bias=='linear'
    pocB(pocB<bias_lim)=(pocB(pocB<bias_lim).*bias_modelB_coef(1)) + bias_modelB_coef(2);
end

% calculate prediction bounds log10(y +/- e), where e is prediction bound defined
% as e = t * sqrt(MSE + xSx'); x is the row vector of the design
% matrix for new predictors ( i.e., [1 log10(bbp) log10(comp)
% log10(comp)*log10(bbp)]; % note that e will be in logspace
x=[ones(size(bbp)) log10(bbp) log10(comp) log10(comp).*log10(bbp)];

for i = 1:length(bbp)
    e(i,1) = t  * sqrt(mse + x(i,:)*S*x(i,:)');
end
% rename final data
POC=pocB;
POC(POC<0)=0; % fix negative values to be 0
pi=[log10(pocB)-e log10(pocB)+e];
PI=10.^pi;
end
