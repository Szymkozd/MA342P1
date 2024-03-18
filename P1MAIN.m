clc
clear variables
close all

numberGantryAngles = 1;
numberPencilBeams = 32; % avoid odd numbers of beams

prescriptionTargetLower = 70;
prescriptionTargetUpper = 95;
prescriptionOAR = 45;
prescriptionTissue = 100;
prescription = [prescriptionTargetLower,prescriptionTargetUpper,prescriptionOAR,prescriptionTissue];

% patientMatrix = patientMatrixModified(15, 0, 20, [-25;0], 20, 1); % patient with tumor and OAR on opposite sides
% patientMatrix = patientMatrixModified(0, 0, 40, [0;0], 30, 1); % patient with OAR inside tumor
% patientMatrix = patientMatrixModified(0, 0, 35, [0;-8], 20, 1); % similar to the example document
% patientMatrix = patientMatrixModified(0, 0, 35, [-30, 0;0, 30], 20, 1); % two OAR on side and top of tumor
% patientMatrix = patientMatrixModified(0, 0, 25, [-30, 0, 30, 0; 0, 30, 0, -30], 20, 1); % four OAR surrounding tumor
patientMatrix = patientMatrixModified(0, 0, 0, [], 0, 1);

optimize(prescription,numberGantryAngles,numberPencilBeams, patientMatrix)
