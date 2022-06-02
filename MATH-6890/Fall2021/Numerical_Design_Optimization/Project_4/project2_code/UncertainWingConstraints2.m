function [c, ceq, dcdx, dceqdx,sts_s] = UncertainWingConstraints2(x, L, E, force, yield, Nelem)
% Computes the nonlinear inequality constraints for the wing-spar problem
% Inputs:
%   x - the DVs; x(1:Nelem+1) inner and x(Nelem+2:2*(Nelem+1) outer radius
%   L - length of the beam
%   E - longitudinal elastic modulus
%   force - force per unit length along the beam axis x
%   yield - the yield stress for the material
%   Nelem - number of finite elements to use
% Outputs:
%   c, ceq - inequality (stress) and equality (empty) constraints
%   dcdx, dceqdx - Jacobians of c and ceq
%--------------------------------------------------------------------------
assert( size(force,1) == (Nelem+1) );
assert( size(x,1) == (2*(Nelem+1)) );

[c,sts_s] = CalcInequality(x);
ceq = [];
dcdx = zeros(2*(Nelem+1),Nelem+1);
dceqdx = [];
% for k = 1:2*(Nelem+1)
%     xc = x;
%     xc(k) = xc(k) + complex(0.0, 1e-30);
%     dcdx(k,:) = imag(CalcInequality(xc))/1e-30;
% end 

    function [cineq,stdDev_stress] = CalcInequality(x)
        % compute the displacements and the stresses
        r_in = x(1:Nelem+1);
        r_out = x(Nelem+2:2*(Nelem+1));
        Iyy = CalcSecondMomentAnnulus(r_in, r_out);
        
        % using a 3 point Gauss-Hermite quadrature
        sigma1 = force(1)/10; sigma2 = force(1)/20;
        sigma3 = force(1)/30; sigma4 = force(1)/40;
        
        xi = [-1.22474487139; 0.0; 1.22474487139];
        wt = [0.295408975151; 1.1816359006; 0.295408975151]./sqrt(pi);
        
        x = (0:L/Nelem:L)';
        
        mean_stress = 0;
        mean_sqstress = 0;
        for i1 = 1:length(xi)
            pt1 = sqrt(2)*sigma1*xi(i1);
            for i2 = 1:length(xi)
                pt2 = sqrt(2)*sigma2*xi(i2);
                for i3 = 1:length(xi)
                    pt3 = sqrt(2)*sigma3*xi(i3);
                    for i4 = 1:length(xi)
                        pt4 = sqrt(2)*sigma4*xi(i4);
                        D = Delta(pt1,pt2,pt3,pt4,L,Nelem);
                        f_u = force + D;
                        u_u = CalcBeamDisplacement(L,E,Iyy,f_u,Nelem);
                        stress_u = CalcBeamStress(L,E,r_out,u_u,Nelem);
                        mean_stress = mean_stress + wt(i1)*wt(i2)*wt(i3)*wt(i4)*stress_u;
                        mean_sqstress = mean_sqstress + wt(i1)*wt(i2)*wt(i3)*wt(i4)*stress_u.*stress_u;
                    end
                end
            end
        end
        %plot(x,f_u); hold on; plot(x,force)
        stdDev_stress = sqrt(mean_sqstress - mean_stress.^2);
        %plot(stdDev_stress)
        %cineq = mean_stress./yield - ones(Nelem+1,1);
        cineq = mean_stress;
    end
end