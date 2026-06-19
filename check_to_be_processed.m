function datatbl=check_to_be_processed
subs=readtable('/data/pt_03187/data/in_vivo/id_info/subject_ids.csv');
subs.bids_subjID=arrayfun(@num2str,subs.bids_subjID,'UniformOutput',false);
for a=1:length(subs.bids_subjID)
    subs.bids_subjID{a}=['sub-',repmat('0',1,3-length(subs.bids_subjID{a})),subs.bids_subjID{a}];
end
%we preload all dirs to loop through them
dirs.bidsdir='/data/pt_03187/data/in_vivo/bids/';
sourcedir='/data/pt_03187/data/in_vivo/source';
dirs.loraks_bids='/data/pt_03187/data/in_vivo/bids/derivatives/LORAKS';
dirs.LCPCA='/data/pt_03187/data/in_vivo/bids/derivatives/LORAKS/derivatives/LCPCA';
dirs.LCPCA_distCorr='/data/pt_03187/data/in_vivo/bids/derivatives/LORAKS/derivatives/LCPCA_distCorr';
dirs.qMRI='/data/pt_03187/data/in_vivo/bids/derivatives/LORAKS/derivatives/qMRI';
dirs.QSMxT='/data/pt_03187/data/in_vivo/bids/derivatives/LORAKS/derivatives/QSMxT';
dirs.nighres='/data/pt_03187/data/in_vivo/bids/derivatives/LORAKS/derivatives/nighres';
dirs.relax_R2='/data/pt_03187/data/in_vivo/bids/derivatives/relax_R2';
dirs.r2prime='/data/pt_03187/data/in_vivo/bids/derivatives/r2prime';
rawdirnames=fieldnames(dirs);
tmpsubs_bids=subs.bids_subjID;

%we grow the final table and since it does not exist the first time we run
%the loop, we have to make a 
firsttry=true;

for a=1:height(subs)
    tmpsub=subs.bids_subjID{a};
    tmptblname=['/data/pt_03187/data/in_vivo/id_info/',subs.bids_subjID{a},'_sessions.csv'];
    
    if exist(tmptblname,'file')~=0
        tmptbl=readtable(tmptblname);
        tmptbl.bids_sesID=arrayfun(@num2str,tmptbl.bids_sesID,'UniformOutput',false);
        tmptbl.sesID=arrayfun(@num2str,tmptbl.sesID,'UniformOutput',false);
        raw_sess_names=mydir(fullfile(sourcedir,subs.subjID{a}),[],2);
        if sum(contains(tmptbl.sesID,raw_sess_names))<length(raw_sess_names)
            %we simply are missing a few raw sessions and will add them
            %later during bidsification. For now, we just give them a
            %running nunber
            add_tbl=cell2table(raw_sess_names(~contains(raw_sess_names,tmptbl.sesID)), ...
                "VariableNames",{'sesID'});
            add_tbl.bids_sesID=arrayfun(@num2str,(1:height(add_tbl))+...
                str2double(tmptbl.bids_sesID{end}),'UniformOutput',false);
            tmptbl=[tmptbl;add_tbl];
        end
        for b=1:length(tmptbl.bids_sesID)
            tmptbl.bids_sesID{b}=['ses-',repmat('0',1,2-length(tmptbl.bids_sesID{b})),tmptbl.bids_sesID{b}];
        end
        tmptbl.bids_subjID=repmat(subs.bids_subjID{a},height(tmptbl),1);
        tmptbl.subjID=repmat(subs.subjID{a},height(tmptbl),1);
    else
        continue
    end
    if firsttry
        tmpsesstbl=tmptbl;
        firsttry=false;
    else
        tmpsesstbl=[tmpsesstbl;tmptbl];
    end
end

rawdirnames=['raw','loraks_processed','dcm','dcm_processed',fieldnames(dirs)'];
datatbl=[tmpsesstbl,cell2table(cell(height(tmpsesstbl),length(rawdirnames)),'VariableNames',rawdirnames)];
datatbl.scanner=cell(height(datatbl),1);
datatbl.scanT=cell(height(datatbl),1);
for a=1:length(rawdirnames)%directories
    tmpdirname=rawdirnames{a};
    for b=1:height(tmpsesstbl)%subject information
        tmpsub=tmpsesstbl.subjID(b,:);
        tmpsub_bids=tmpsesstbl.bids_subjID(b,:);
        tmpsess=tmpsesstbl.sesID{b};
        tmpsess_bids=tmpsesstbl.bids_sesID{b};
        switch tmpdirname
            case 'raw'
                doi=fullfile(sourcedir,tmpsub,tmpsess,'raw');
            case 'loraks_processed'
                doi=fullfile(sourcedir,tmpsub,tmpsess,'nii_loraks_recon');
            case 'dcm'
                doi=fullfile(sourcedir,tmpsub,tmpsess,'dcm');
            case 'dcm_processed'
                doi=fullfile(sourcedir,tmpsub,tmpsess,'nii');
            case 'bidsdir'
                doi=fullfile(dirs.bidsdir,tmpsub_bids, ...
                   tmpsess_bids,'anat');
            case 'loraks_bids'
                doi=fullfile(dirs.bidsdir,'derivatives','LORAKS',tmpsub_bids, ...
                    tmpsess_bids,'anat');
            otherwise
                doi=fullfile(dirs.bidsdir,'derivatives','LORAKS','derivatives', ...
                    tmpdirname,tmpsub_bids, ...
                    tmpsess_bids,'anat');

        end
        if ~exist(doi,'dir')&isempty(mydir(doi))
            if         strcmp(tmpdirname,'dcm_processed')
                doi=fullfile(sourcedir,tmpsess,'dcm2niix');
            elseif strcmp(tmpdirname,'QSMxT')
                doi=mydir(fullfile(dirs.bidsdir,'derivatives','LORAKS','derivatives', ...
                    tmpdirname,'/*/',tmpsub_bids, ...
                    tmpsess_bids,'anat'),[],1);
                if ~isempty(doi)
                doi=fileparts(doi{1});
                end
            elseif contains(lower(tmpdirname),'r2')
                doi=fullfile(dirs.bidsdir,'derivatives', ...
                    tmpdirname,tmpsub_bids, ...
                    tmpsess_bids,'anat');
            end
           
        end
        if  exist(doi,'dir')~=0&~isempty(mydir(doi))
            datatbl.(tmpdirname){b}=doi;
            if strcmp(tmpdirname,'loraks_processed')
                data=mydir([doi,'/*mtw_*.json'],[],1);
                if isempty(data)
                      datatbl.(tmpdirname){b}=[];
                    continue
                end
                json=readstruct(data{1});
                datatbl.scanner{b}=json.ManufacturersModelName;
                json.ManufacturersModelName;
                if strcmp(datatbl.scanner{b},'MAGNETOM Terra.X')
                    datatbl.scanT{b}='7';
                elseif strcmp(datatbl.scanner{b},'Prisma_fit')
                    datatbl.scanT{b}='3';
                 end
            end
        end
       

    end
end
for a=1:height(datatbl)
    if ~isempty(datatbl.scanner{a})
        if contains(datatbl.scanner{a},'magnetom','IgnoreCase',true)
            datatbl.scanT{a}='7';
        else
            datatbl.scanT{a}='3';
        end
    end
end
end