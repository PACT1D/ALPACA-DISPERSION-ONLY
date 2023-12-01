% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% Sparse Data Definition File                                      
%                                                                  
% Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  
%       (http://www.cs.vt.edu/~asandu/Software/KPP)                
% KPP is distributed under GPL, the general public licence         
%       (http://www.gnu.org/copyleft/gpl.html)                     
% (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           
% (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            
%     With important contributions from:                           
%        M. Damian, Villanova University, USA                      
%        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
%                                                                  
% File                 : mech_Sparse.m                             
% Time                 : Tue Jan 14 15:47:22 2020                  
% Working directory    : /home/katie/PACT1D/PACT-1D-CalNex_KBTworking/mechanism
% Equation file        : mech.kpp                                  
% Output root filename : mech                                      
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




%  ----------> Sparse Jacobian Data                                

% LU_IROW - Row indexes of the LU Jacobian of variables            
 global LU_IROW;
% LU_ICOL - Column indexes of the LU Jacobian of variables         
 global LU_ICOL;
% LU_CROW - Compressed row indexes of the LU Jacobian of variables 
 global LU_CROW;
% LU_DIAG - Diagonal indexes of the LU Jacobian of variables       
 global LU_DIAG;

