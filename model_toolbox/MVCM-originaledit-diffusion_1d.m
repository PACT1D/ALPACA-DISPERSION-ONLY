%1D diffusion implemented by Jennie Thomas & Cyril Falvo
%Edited by Jochen Stutz - 18 Sept 2019
%Edited by Meeta Cesler-Maloney October 2022 to add mixing with a background column after vertical diffusuion
%Edited by Meeta Cesler-Maloney March 2023 fixing vertical and horizontal transport rate calculations

% MVCM 03/23: added new variables for horizontal transport rates and concentration of species exported by horizontal advection
function [ spec_t1 , VT_t , HT_t , depo_t, total_loss_to_ground] = diffusion_1d (spec_t0, spec_bkg, BOX_WALL, BOXCH, ...
    Kz_interp, exch_coeff_interp, diff_const_interp, rho, temperature, NLEV, ...,
    dt_diff_in, dt_diff_horiz, Eff_dep_vel_t,  n_step_diff, total_loss_to_ground, run_chem, add_surface_source_HONO, horz_mix)

%Do 1 dimensional diffusion for this time step, takes in the current
%concentration of species at t=0 and returns species concentrations at t=t+1

mech_Parameters;

% MVCM 03/23: added arrays to store seperated vertical and horizontal transport rates and concentrations of species exported by horizontal advection
%set up array to save before and after diffusion concentrations
% species concentration before any transport
spec_before_diff       = zeros(NVAR,NLEV);
% species concentration after vertical transport only
spec_after_Vt        = zeros(NVAR,NLEV);
% species concentration after horizontal transport only
spec_after_Ht       = zeros(NVAR,NLEV);
% species concentration exported during horizontal advection

%set up array to save deposition rate, vertical transport rate and surface sources
% MVCM added array to store horizontal transport rates
VT_t                       = zeros(NVAR,NLEV); % vertical transport rate
HT_t                       = zeros(NVAR,NLEV); % horizontal transport rate
depo_t                     = zeros(NVAR,1);     %the total deposition rate
depo_t_inst                = zeros(NVAR,1);     %sub_timestep deposition rate, for surface source
loss_to_ground_inst        = zeros(NVAR,1);     %sub timestep loss to ground for each species in molecules/cm2
diff_spec                  = zeros(NVAR,NLEV);  %difference during diffusion for each species due to diffusion/deposition on model levels
sum_diff                   = zeros(NVAR,1);     %vertical column change before and after diffusion/deposition

%save box wall to wall distance in cm
BOX_WALL_cm = BOX_WALL*100. ; %box wall to wall distance in cm

%how much backward vs. forward integraton to consider
%backward alpha = 1.0, forward  alpha = 0.0
alpha = 1.0;

%set the diffusion time step should be set according to the lowest box
%Kz and the time step such that the following is true
%N_TIME_DIFF = ceil(dt_diff_in*Kz(1)/BOX_WALL(1)/(BOXCH(2)-BOXCH(1)));
%N_TIME_DIFF = 1;

%disp(['running diffusion for ', num2str(N_TIME_DIFF) ,' steps']);

%dt_diff = dt_diff_in/N_TIME_DIFF;
dt_diff = dt_diff_in;
%set spec t1 = spec t0
spec_t1 = spec_t0;

%Calculate online deposition velocity if run chem is on
if ((run_chem == 1) && (add_surface_source_HONO == 1))
    %Calcualte online deposition velocity - modify online_depo_vel
    % to change deposition velocities during a run
    Eff_dep_vel_t = online_depo_vel(Eff_dep_vel_t);
end

%loop over diffusion time steps
%for ii = 1:N_TIME_DIFF
    %loop over species
    for nn = 1:NVAR
        
        %Get the effective deposition velocity for this species
        v_eff_spec = Eff_dep_vel_t(nn);
        
        %Interpolate the Kz values in time
        %contrib_t1 = (N_TIME_DIFF-(ii-1))/N_TIME_DIFF;
        %Kz = Kz1*contrib_t1+Kz2*(1.-contrib_t1);
        
        %Interpolate Diffusion Constant in Time
        %diffusion_constant=diff_const1(nn,:)*contrib_t1+diff_const2(nn,:)*(1.-contrib_t1);       
                
        %Calcuate total diffusion constant - Kz + molecular diffusion
        Kz = Kz_interp + diff_const_interp';
        
        %DO 1D diffusion and deposition
        %save concentration prior to all diffusion
        spec_before_diff(nn,:) = spec_t1(nn,:);
        %convert to mixing ratio
        for k=1:NLEV
            spec_t1(nn,k) = spec_t1(nn,k)/rho(k); 
        end
        %calcuate the first coefficients needed
        [a_coeffs, b_coeffs, c_coeffs, dp, dm, beta] = coeffs_diffusion_1d_first_step(Kz, rho, BOX_WALL, BOXCH, alpha, dt_diff, NLEV, v_eff_spec);
        %calcuate the remaining coeffiencts
        [d_coeffs] = coeffs_diffusion_1d(spec_t1(nn,:),beta,dp,dm,NLEV, v_eff_spec, dt_diff,BOX_WALL);
        %do the tridiag algorithim
        spec_t1(nn,:) = solve_triadiag(a_coeffs, b_coeffs, c_coeffs, d_coeffs, NLEV);
        % save species concentration after vertical transport
        spec_after_Vt(nn,:) = spec_t1(nn,:);
        % convert species concentration after vertical transport only back to concentration
        for k=1:NLEV
            spec_after_Vt(nn,k) = spec_t1(nn,k) * rho(k); 
        end
        % MVCM 05/2023: uncomment if statement to manually check timing
        %if nn == 1
        %    disp('doing vertical diffusion');
        %end
        %
        % MVCM 03/23: added horizontal advection process below
        % for each species, acccount for horizontal advection and convert back to number concentration
        for k=1:NLEV
            if horz_mix == 1
                % MVCM 05/2023: uncomment if statement to manually check timing
                %if k == 1 && nn == 1
                %    disp('doing horizontal diffusion');
                %end
                %
                % multiply exchange coeff (s-1) * horizontal mixing timestep (s) to get unitless fraction background and fraction polluted
                frac_bkg_step = exch_coeff_interp(k) * double(dt_diff_horiz);
                frac_polluted = 1 - frac_bkg_step;
                % calculate species mixing ratio  after horizontal advection
                spec_t1(nn,k) = (frac_bkg_step * spec_bkg(nn,k)) + (frac_polluted * spec_t1(nn,k));
                % convert all species mixing ratios back to concentrations
                spec_after_Ht(nn,k) = spec_t1(nn,k) * rho(k);
                spec_t1(nn,k) = spec_t1(nn,k) * rho(k);
            else
                % if horizontal advection is not happening, leave spec_t1 as is and fill others with zeroes
                spec_after_Ht(nn,k) = double(0);
                spec_t1(nn,k) = spec_t1(nn,k) * rho(k);
            end
        end
        %END 1D diffusion, horizontal advection and deposition
    end %end loops over species


    %calculate deposition rate for this time
    %in units of molec/cm2/s
    % MVCM 03/23: uses the species after vertical transport only, consider changing this later
    for nn=1:NVAR
        diff_spec(nn,1)=spec_before_diff(nn,1)*BOX_WALL_cm(1) - spec_after_Vt(nn,1)*BOX_WALL_cm(1);
        for k=2:NLEV
            diff_spec(nn,k)=spec_before_diff(nn,k)*(BOX_WALL_cm(k)-BOX_WALL_cm(k-1)) - spec_after_Vt(nn,k)*(BOX_WALL_cm(k)-BOX_WALL_cm(k-1));
        end
        %save the deposition rate for this instant
        
        sum_diff(nn)=sum(diff_spec(nn,:));
        %save and sum the total deposition rate, divide by the timestep at the very end
        loss_to_ground_inst(nn) = sum_diff(nn);     %positive term, ground storage
        depo_t(nn)            = depo_t(nn) + sum_diff(nn);  %depo_t is the total deposition, summed here - divided by time step later
    end
    
    %make cumulative surface storage sum
    total_loss_to_ground = total_loss_to_ground+loss_to_ground_inst  ; %molec/cm2
    
%end  %end loop over sub timesteps

%calculate vertical trasnport rate - this is change in number/cm3/s for each level
%This is due to vertical mixing and deposition
%VT_t = (spec_t1-spec_t0)/double(dt_diff_in);

% MVCM 03/23: updated vertical and horizontal transport rate calculations
VT_t = (spec_after_Vt-spec_t0)/double(dt_diff_in);

if horz_mix == 1
    HT_t = (spec_after_Ht-spec_after_Vt)/double(dt_diff_horiz);
else
    HT_t = double(0);
end


%Divide deposition amount by the timestep, number/cm2/s
depo_t = depo_t/double(dt_diff_in);

end
