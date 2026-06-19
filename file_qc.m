function filelist = file_qc(filelist)
dirs=filelist.Properties.VariableNames;
dirs(contains(dirs,'ses')|contains(dirs,'sub')|contains(dirs,'scan'))=[];
for a=1:length(dirs)
    tmpdir=dirs{a};
    for b=1:height(filelist)
        tmpsub=filelist.bids_subjID(b,:);
        tmpses=filelist.bids_sesID{b};
        basepath=cell2mat(filelist{b,tmpdir});
        if isempty(filelist{b,tmpdir})
            continue
        end
        switch tmpdir
            case 'raw'
            case 'loraks_processed'
                %do we have all files, are they not noisy?
            case 'dcm'
                %do we have all files, are they properly processed?
            case 'dcm_processed'
                
            case 'bidsdir'
                %do we have all files, are they properly bidsified?
            case 'loraks_bids'
                %do we have all files, are they properly bidsified?
            case 'LCPCA'
            case 'LCPCA_distCorr'
            case 'qMRI'                
                %Did we take the correct input files?
                jsonstring=jsondecode(fileread(fullfile(basepath,[tmpsub,'_',tmpses,'_PDmap.json'])));
                qmriorgdir=fileparts(jsonstring.history.input(1).filename);
                if contains(qmriorgdir,tmpses)&&contains(qmriorgdir,tmpsub)&&contains(qmriorgdir,'distCorr')
                else
                    
                end
            case 'QSMxT'
            case 'nighres'
            case 'relax_R2'
            case 'r2prime'
            otherwise 
                error(['The temporary directory ',tmpdir,' was not recognized. Check it again.'])
        end
    end
end
end

