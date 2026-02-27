% This function file generates the xdmf file needed to visualize the h5 file in paraview 

% function XDMF_TEXT_GENERATOR(filename,h5string,h5att_name)

function XDMF_TEXT_GENERATOR(filename,h5string,h5att_name,z_depth,y_breadth,x_breadth)

fid = fopen(filename,'w+t');
if fid <0;
    fprintf('error openingfile\n');
    return
end

% z_depth = 2560;
% x_breadth = 2560;
% y_breadth = 299;

fprintf(fid,'%s\n','<?xml version="1.0"?>');
fprintf(fid,'%s\n','<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>');
fprintf(fid,'%s\n','<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">');
fprintf(fid,' %s\n','<Domain>');
fprintf(fid,'\n');
fprintf(fid,'  %s\n','<Grid Name="Cell Data" GridType="Uniform">');
fprintf(fid,'    %s%d %d %d%s\n','<Topology TopologyType="3DCoRectMesh" Dimensions="',z_depth+1, x_breadth+1, y_breadth+1,'"></Topology>');
fprintf(fid,'    %s\n','<Geometry Type="ORIGIN_DXDYDZ">');
fprintf(fid,'      %s\n','<!-- Origin -->');
fprintf(fid,'      %s\n','<DataItem Format="XML" Dimensions="3">0 0 0</DataItem>');
fprintf(fid,'      %s\n','<!-- DxDyDz (Spacing/Resolution)-->');
fprintf(fid,'      %s\n','<DataItem Format="XML" Dimensions="3">1 1 1</DataItem>');
fprintf(fid,'    %s\n','</Geometry>');
fprintf(fid,'\n');

% fprintf(fid,'      %s\n','<Attribute Name="UNDEFORMED" AttributeType="Scalar" Center="Cell">');
% fprintf(fid,'      %s%d %d %d%s\n','<DataItem Format="HDF" Dimensions="',z_depth, x_breadth, y_breadth,'" NumberType="uint8" Precision="4"');
% % h5string = sprintf('%s','FIB_.h5:/Data/DEF');
% fprintf(fid,'       %s\n',h5string);
% fprintf(fid,'      %s\n','</DataItem>');
% fprintf(fid,'       %s\n','</Attribute>');

for ii = 1:length(h5string)
fprintf(fid,'      %s%s%s\n','<Attribute Name="',h5att_name{ii},'" AttributeType="Scalar" Center="Cell">');
fprintf(fid,'      %s%d %d %d%s\n','<DataItem Format="HDF" Dimensions="',z_depth, x_breadth, y_breadth,'" NumberType="double" Precision="4" >');
% h5string = sprintf('%s','FIB_.h5:/Data/DEF');
fprintf(fid,'       %s\n',h5string{ii});
fprintf(fid,'      %s\n','</DataItem>');
fprintf(fid,'       %s\n','</Attribute>');
if ii < length(h5string)
fprintf(fid,'\n');
end
end

fprintf(fid,'\n');
fprintf(fid,'  %s\n','</Grid>');
fprintf(fid,'    %s\n','<!-- *************** END OF Cell Data *************** -->');
fprintf(fid,'\n');
fprintf(fid,' %s\n','</Domain>');
fprintf(fid,'%s\n','</Xdmf>');

fclose(fid);







