%% Calculate new fluxes for epimetoncofit
% @author: Scott Campit

initCobraToolbox;
changeCobraSolver('gurobi');

% Load model 
model_path = 'C:\Users\scampit\Desktop\epimetoncofit\gssm\';
model = readCbModel([model_path, 'A498.xml'],'fileType', 'SBML');

