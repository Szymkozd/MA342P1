function plotDMatrix(D,patientMatrix,x, plotContour,plotMesh, titleText)
D = sum(D*x,2);
[row,col] = find(patientMatrix'>0);
Dose = zeros(size(patientMatrix,1),size(patientMatrix,2));
for i = 1:length(row)
    Dose(row(i),col(i)) = D(i);
end
if(plotContour)
    figure()
    set(gcf, 'Position', [100, 50, 800, 600]);
    mesh(Dose)
    title(titleText)
    axis([200,300, 200,300])
end
if(plotMesh)
    figure()
    hold on
    imagesc(Dose);
    contour(patientMatrix, 'LineColor', 'w', 'LineWidth', 3);
    axis([200,300, 200,300])
    title(titleText)
    hold off
end
end
