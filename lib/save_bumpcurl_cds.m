function [cds] = save_bumpcurl_cds(input_struct)

% Make CDS file
cds = commonDataStructure();
cds.file2cds([input_struct.folder filesep input_struct.fname],input_struct.ranBy,input_struct.array,input_struct.monkey,input_struct.lab,'ignoreJumps',input_struct.task,input_struct.mapfile);

if(isfield(input_struct,'save_folder'))
    save(input_struct.save_folder,'cds','-v7.3');
end