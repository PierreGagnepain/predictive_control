function [pos_vx] = mm2vx_mat(pos_mm,M);
%------------------------------------------------------------------------
% FORMAT function [pos_vx] = mm2vx_mat(pos_mm,M);
% Transforms point(s) 'pt_mm' mm coordinates into voxel coordinates 'pt_vx'
% according to the transformation matrix 'M'.
% 'M' is the 4x4 affine transformation matrix : from vx to mm
%------------------------------------------------------------------------
	
	if size(pos_mm,1)==3
        	Npt = size(pos_mm,2);
	elseif size(pos_mm,2)==3;
        	pos_mm = pos_mm';
        	Npt = size(pos_mm,2);
	else
        	error('Wrong vectors format !')
	end
	
	tmp = M\[pos_mm;ones(1,Npt)];
	pos_vx = tmp(1:3,:);
