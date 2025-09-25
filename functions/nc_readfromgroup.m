function data = nc_readfromgroup(fname, group, varname, varargin)
% READFROMGROUP  Read a NetCDF4 variable from a specific group.
%   data = readFromGroup(fname, group, varname)
%   data = readFromGroup(fname, group, varname, start, count)
%
% Example:
%   HSo = readFromGroup('shorelines_output.nc', 'projected_grid_003','HSo');
  
  ncid  = netcdf.open(fname,'NC_NOWRITE');
  cleanupObj = onCleanup(@() netcdf.close(ncid));
  gid   = netcdf.inqNcid(ncid, group);
  vid   = netcdf.inqVarID(gid, varname);
  if nargin==3
    data = netcdf.getVar(gid, vid);
  else
    start = varargin{1};
    count = varargin{2};
    data  = netcdf.getVar(gid, vid, start, count);
  end
end
