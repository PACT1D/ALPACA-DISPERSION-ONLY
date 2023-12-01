% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% Auxiliary Routines File                                          
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
% File                 : mech_Util.m                               
% Time                 : Tue Jan 14 15:47:22 2020                  
% Working directory    : /home/katie/PACT1D/PACT-1D-CalNex_KBTworking/mechanism
% Equation file        : mech.kpp                                  
% Output root filename : mech                                      
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




% User INLINED Utility Functions                                   

% End INLINED Utility Functions                                    

% Utility Functions from KPP_HOME/util/util                        
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% UTIL - Utility functions                                         
%   Arguments :                                                    
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ****************************************************************
%                            
% InitSaveData - Opens the data file for writing
%
% ****************************************************************

function InitSaveData ()

global mech_FID

      mech_FID = fopen('mech.dat','w');

return %  InitSaveData

% End of InitSaveData function
% ****************************************************************

% ****************************************************************
%                            
% SaveData - Write LOOKAT species in the data file 
%
% ****************************************************************

function SaveData ()

global VAR FIX CFACTOR LOOKAT NLOOKAT mech_FID

      C(1:142) = VAR(1:142);
      C(142+1:148) = FIX(1:6);
      
      fprintf(mech_FID,'%12.5e,',C(LOOKAT(1:NLOOKAT)));

return %  SaveData

% End of SaveData function
% ****************************************************************

% ****************************************************************
%                            
% CloseSaveData - Close the data file 
%
% ****************************************************************

function CloseSaveData ()
global mech_FID

      fclose( mech_FID );

return %  CloseSaveData

% End of CloseSaveData function
% ****************************************************************


% End Utility Functions from KPP_HOME/util/util                    
% End of UTIL function                                             
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


