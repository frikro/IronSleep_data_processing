there=@(x)~cellfun(@isempty,x);
filelist=check_to_be_processed;
%filelist=file_qc(filelist);
addpath('/data/u_krohn_software/git/bidsification_mpm')
filelist(~there(filelist.loraks_processed),:)=[];
missing_bidsification=filelist(~there(filelist.bidsdir),:);
sourcedir='/data/pt_03187/data/in_vivo/source';
% for a=1:height(missing_bidsification)
%     dcm2nii_matlab_version(missing_bidsification.dcm{a})
% end
% for a=1:height(filelist)
%     if ~exist(fullfile(fileparts(filelist.loraks_bids{a}),'fmap'),'dir')
%         system(['ln -s ',fullfile(fileparts(filelist.bidsdir{a}),'fmap'),' ',...
%             fullfile(fileparts(filelist.loraks_bids{a}),'fmap')])
%     end
% end
no_nighres=filelist(~there(filelist.nighres)&there(filelist.loraks_bids)&~there(filelist.nighres),:);
no_qmri=filelist(~there(filelist.qMRI)&there(filelist.loraks_bids)&~there(filelist.nighres),:);
%writetable(filelist,'overview_table.csv')
%writetable(filelist,'overview_table.xlsx')
mass_distributed_process_generalized_histopark(no_qmri.scanT,'loraks','complex','/data/pt_03187/data/in_vivo/bids/derivatives/LORAKS/',no_nighres.bids_subjID,no_nighres.bids_sesID)