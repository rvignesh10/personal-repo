%% Calculate the mean statistics of stresses and displacement due to uncertainties

function [mStress,stdDevStress,m_u,stdDevu] = uncertainity(zmax,f_nom,Iyy,E,L,Nelem)
% Inputs - zmax   - R_out values at each node point - array of size [Nelem+1 ,1]
%          f_nom  - nominal force values at the node locations
%          Iyy    - Area Moment of Inertia at each node location
%          E      - Young's Modulus
%          L      - Length of the Spar
%          Nelem  - Num of elements 
% Outputs- mStress- mean normal stress due to uncertain loading at each
%                   nodal location      
%     stdDevStress- standard deviation of normal stress due to uncertain
%                   loading at nodal location
%             m_u - Mean Displacement due to uncertain loading
%          stdDevu- standard deviaiton of displacement due to uncertain
%                   loading

% using a 3 point Gaussian quadrature rule . 
xi = [-1.22474487139; 0.0; 1.22474487139];
wt = [0.295408975151; 1.1816359006; 0.295408975151]./sqrt(pi);
% standard deviation of the perturbation variables . 
sigma1 = f_nom(1)/10; sigma2 = f_nom(1)/20;
sigma3 = f_nom(1)/30; sigma4 = f_nom(1)/40;

mStress = 0;
m_sqStress = 0;
m_u = 0;
m_u2 = 0;
for i1 = 1:length(xi)
    pt1 = sqrt(2)*sigma1*xi(i1);
    for i2 = 1:length(xi)
        pt2 = sqrt(2)*sigma2*xi(i2);
        for i3 = 1:length(xi)
            pt3 = sqrt(2)*sigma3*xi(i3);
            for i4 = 1:length(xi)
                pt4 = sqrt(2)*sigma4*xi(i4);
                D = Delta(pt1,pt2,pt3,pt4,L,Nelem);
                f_u = f_nom + D;
                u_u = CalcBeamDisplacement(L,E,Iyy,f_u,Nelem);
                stress_u = CalcBeamStress(L,E,zmax,u_u,Nelem);
                % compute mean stress
                mStress = mStress + wt(i1)*wt(i2)*wt(i3)*wt(i4)*stress_u;
                % compute mean stress square
                m_sqStress = m_sqStress + wt(i1)*wt(i2)*wt(i3)*wt(i4).*stress_u.*stress_u;
                % compute mean displacement
                m_u = m_u + wt(i1)*wt(i2)*wt(i3)*wt(i4)*u_u;
                % compute mean displacement square
                m_u2 = m_u2 + wt(i1)*wt(i2)*wt(i3)*wt(i4)*u_u.*u_u;
            end
        end
    end
end
% calculate standard deviaiton of stress
stdDevStress = sqrt(m_sqStress - mStress.^2);

% calculate standard deviation of displacement
stdDevu = sqrt(m_u2 - m_u.^2);
end