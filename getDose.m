function D = getDose(patientMatrix, numberGantryAngles,numberPencilBeams)
    offAxisData = readmatrix('OffAxisData.txt');
    expSampW = offAxisData(:,1);
    expSampO = offAxisData(:,2);
        
    [x,y] = find(patientMatrix>0);
    
    D = zeros(size(x,1),numberGantryAngles*numberPencilBeams);
    
    for alpha = 1:numberGantryAngles
        a = ((2*pi)/numberGantryAngles)*(alpha-1);
        for i = 1:numberPencilBeams
            for j = 1:size(x,1)
                p = [x(j)-250,y(j)-250]*0.2; %recenter at isocenter
                calcPenMod = pencilModel(p,a,i,expSampW,expSampO,numberPencilBeams);
                D(j,((alpha-1)*numberPencilBeams)+i) = calcPenMod;
            end
        end
    end
end

% save('DoseMatrix.mat',"D")

function D = pencilModel(p,a,i,expSampW,expSampO,numPen)
    P_0 = 0.873;
    mu = 0.0469;
    r = 1.683;
    gamma = 5.2586;
        
    l_gc = 100; % distance from the gantry to the isocenter (cm)

    [o,d,s] = pencilDistCalc(p,a,i,numPen);
    
    alpha_d = 0.1299-0.0306*log(d);
    l_gp = sqrt((s + d)^2 + o^2);
    ISF = (l_gc/l_gp)^2;
    psi = (abs(o)*l_gc)/(s + d);
    O = offAxisFactor(psi, expSampW, expSampO);

    if d >= 1.5
        D = (P_0*exp(-mu*(d-1.5))*(1-exp(-gamma*r))+((r*d*alpha_d)/(r+1.5)))*ISF*O;
    elseif (0 <= d && d < 1.5)
        D = (0.4*(d/1.5)+0.6)*(P_0*(1-exp(-gamma*r))+((1.5*r*alpha_d)/(r+1.5)))*ISF*O;
    else
        D = 0;
    end
end

function [o,d,s] = pencilDistCalc(p,a,i,numPen)
    l_gc = 100; %cm
    rot_mat = [cos(a), -sin(a); sin(a), cos(a)];
    cbp = rot_mat*[0; -((i*0.3)-((numPen/2)*0.3)-0.15)];
    cbp_vec = [-l_gc; -((i*0.3)-((numPen/2)*0.3)-0.15)]; % for even pencils, up most is i=1
    cbp_vec = rot_mat*cbp_vec;
    g = [l_gc*cos(a), l_gc*sin(a)]; %need to recenter to at gantry

    gp_vec = p - g;

    sd = (dot(gp_vec,cbp_vec)/(norm(cbp_vec))); %scalar projection
    o = real(sqrt(norm(gp_vec)^2 - sd^2));

    a = (g(1) - cbp(1))^2 + (g(2) - cbp(2))^2;
    b = 2*(g(1) - cbp(1))*(cbp(1)) + 2*(g(2) - cbp(2))*(cbp(2));
    c = (cbp(1))^2 + (cbp(2))^2 - 10^2;

    t = (2*c)/(-b+sqrt(b^2-4*a*c));

    intsecx = cbp_vec(1)*t;
    intsecy = cbp_vec(2)*t;

    s = norm([g(1)-intsecx,g(2)-intsecy(1)]);

    d = sd - s;
end

function O = offAxisFactor(psi, expSampW, expSampO)
    if psi > 2.5
        O = 0;
    else
        O = interp1(expSampW, expSampO, psi);
    end
end