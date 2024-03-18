function DVH(Dose, x, targetDoseRowIndex, prescriptionTargetLower, titleText)
    Dose = Dose*x./prescriptionTargetLower;
    D_T = Dose(targetDoseRowIndex==2,:);
    D_C = Dose(targetDoseRowIndex==3,:);
    D_N = Dose(targetDoseRowIndex==1,:);
    percentTargetDose = linspace(0, 1.40, 100);
    for i = 1:length(percentTargetDose)
        DVH_T(i) = sum(D_T>percentTargetDose(i));
        DVH_C(i) = sum(D_C>percentTargetDose(i));
        DVH_N(i) = sum(D_N>percentTargetDose(i));
    end
    DVH_T = DVH_T./size(D_T, 1);
    DVH_C = DVH_C./size(D_C, 1);
    DVH_N = DVH_N./size(D_N, 1);
    figure()
    hold on
    linewidth = 3;
    plot(percentTargetDose,DVH_T, 'r-', 'LineWidth', linewidth);
    plot(percentTargetDose,DVH_C, 'g-', 'LineWidth', linewidth);
    plot(percentTargetDose,DVH_N, 'b-', 'LineWidth', linewidth);
    xlabel("Percent Target Dose")
    ylabel("Percent Volume")
    xline(1);
    title(titleText);
    hold off
end