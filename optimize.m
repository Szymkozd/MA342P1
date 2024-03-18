function optimize(prescription,numberGantryAngles,numberPencilBeams, patientMatrix)
% load('patientMatrix.mat');
% load('DoseMatrix4.mat');

prescriptionTargetLower = prescription(1);
prescriptionTargetUpper = prescription(2);
prescriptionOAR = prescription(3);
prescriptionTissue = prescription(4);

D = getDose(patientMatrix, numberGantryAngles,numberPencilBeams);

plotDMatrix(D,patientMatrix,ones(size(D,2)), 1, 1, "Initial dose per fluence for " + numberGantryAngles + " gantry angles and " + numberPencilBeams + " pencil beams");

numberTargetPoints = sum(patientMatrix(:)==2);
numberOARPoints = sum(patientMatrix(:)==3);
numberTissuePoints = sum(patientMatrix(:)==1);
totalPatientPoints = numberTargetPoints + numberOARPoints + numberTissuePoints;

e_Target = ones(numberTargetPoints,1);
e_OAR = ones(numberOARPoints,1);
e_Tissue = ones(numberTissuePoints,1);

%%
% rows of dose matrix corresponding to each type of tissue
targetDoseRowIndex = nonzeros(patientMatrix); 

D_T = D(targetDoseRowIndex==2,:);
D_C = D(targetDoseRowIndex==3,:);
D_N = D(targetDoseRowIndex==1,:);

% upper and lower prescriptions
L_T = prescriptionTargetLower*e_Target;
U_T = prescriptionTargetUpper*e_Target;
U_C = prescriptionOAR*e_OAR;
U_N = prescriptionTissue*e_Tissue;

nu = 25/45;

lambda_T = 2000/numberTargetPoints;
lambda_C = 1800/numberOARPoints;
lambda_N = 100/numberTissuePoints;

zRows = numberPencilBeams*numberGantryAngles + numberTargetPoints + numberOARPoints + numberTissuePoints;

A = [-D_T, zeros(size(D_T,1),  (zRows - size(D_T,2)));
      D_T, -eye(size(D_T,1)),  zeros(size(D_T,1),zRows - size(D_T,1) - size(D_T,2));
      D_C, zeros(size(D_C,1),size(D_T,1)), -eye(size(D_C,1)), zeros(size(D_C,1),zRows - size(D_C,2) - size(D_T,1) - size(D_C,1));
      D_N, zeros(size(D_N,1), zRows-size(D_N,2)-size(D_N,1)), -eye(size(D_N,1))];

b = [-L_T;
      U_T;
      U_C;
      U_N];

c = [zeros(zRows-totalPatientPoints,1); % was totalPatientPoints
    lambda_T*e_Target;
    lambda_C*e_OAR;
    lambda_N*e_Tissue];

lb = [zeros(1,numberPencilBeams*numberGantryAngles), zeros(1,numberTargetPoints), -nu*U_C', zeros(1,numberTissuePoints)]';

[z, optVal] = linprog(c, A, b, [], [], lb, []);
x = z(1:numberPencilBeams*numberGantryAngles);

plotDMatrix(D,patientMatrix,x, 1, 1, "Optimized dose for " + numberGantryAngles + " gantry angles and " + numberPencilBeams + " pencil beams")
DVH(D,x, targetDoseRowIndex, prescriptionTargetLower, "DVH for " + numberGantryAngles + " gantry angles and " + numberPencilBeams + " pencil beams")

% z = [x', gamma', rho', sigma']';
% x is fluence map, gamma is the overdose to the target, sigma is overdose
% to the normal tissues, rho is overdose to the OARs
end
