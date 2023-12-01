% script added by MVCM 03/23
% code copied from other read input file scripts 
function [bkg_in,bkg_list_input,error_msg] = read_bkg(bkgfile,BOX_WALL,Times);
% read background column concentrations in ppb

%set error message to nothing
error_msg = [];

%add model parameters
mech_Parameters;

%-------------------------  Read background file -------------------------------------
%read model levels from the background file
ncid = netcdf.open(bkgfile,'NOWRITE');

 %get the times of the initial concentration file
 varid      = netcdf.inqVarID(ncid,'Times');
 Times_bkg = netcdf.getVar(ncid,varid);
 Times_bkg = Times_bkg.';
 
 %there is on more time in the master time file, compare with the size
 %ntimes-1
 [ntimes ~] = size(Times);

 %if the initial concentraiton time is not equal to the master time file
 %initial time, exit with an error message
 if Times(1:ntimes,:) ~= Times_bkg
   error_msg = ['Times bkg_col.nc are not the same as in the master times file master_time_lev.nc.'];
   return
 end
 
 %get the Box hights from the initial concentration file
 varid      = netcdf.inqVarID(ncid,'BOX_WALL');
 BOX_WALL_bkg = netcdf.getVar(ncid,varid);
 
 %if the model levels are different - then  send an error message
 if BOX_WALL_bkg ~= BOX_WALL 
   error_msg = ['The model levels in the master initialization file master_time_lev.nc are different from the J_values.nc file.'];
   return
 end
 
 %get the names of the bkg in the files
 %get list of varIDs in the files
 varIDs = netcdf.inqVarIDs(ncid);
 [~, nvarIDs] = size(varIDs);
 
 [NLEV ~] = size(BOX_WALL);

 %make space for bkg
 bkg_in = zeros(nvarIDs-1,NLEV,ntimes);  %-1 here, arrays start at 1, there are two varids not used BOX_WALL and Times
 %bkg_in = zeros(nvarIDs,NLEV,ntimes);

 %save bkg and a list of the bkg names
 bkg_list_tmp = {};
for i=0:nvarIDs-1  %Loop over all variables in netcdf file
%for i=0:nvarIDs
  var = netcdf.inqVar(ncid,i);                              %species name from the netcdf file
  if ismember(var,'BOX_WALL')
        disp('reading model levels');
    elseif ismember(var,'Times')
        disp('reading time values');
  else
    bkg_list_tmp = [bkg_list_tmp{:} {var}];
    disp(['reading background species - ' var]);   %finding this species name in the mechanism
    %size(netcdf.getVar(ncid,varIDs(i)))
    %size(bkg_in(i-1,:,:))
    %bkg_in(i-1,:,:) = netcdf.getVar(ncid,i);     
    temp = netcdf.getVar(ncid,i);
    %bkg_in(i-1,:,:) = netcdf.getVar(ncid,i)';  %JPS why do we now have to transpose the input matrix
    bkg_in(i-1,:,:) = netcdf.getVar(ncid,i);
  end
end

 bkg_list_input = bkg_list_tmp;

 return
