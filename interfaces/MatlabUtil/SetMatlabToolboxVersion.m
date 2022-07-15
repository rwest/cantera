function SetMatlabToolboxVersion(ver)
    % Sets the version of the Cantera Matlab Interface
    %
    % :param ver:
    %    string  'legacy'/'new';
    
    % Revert search path
    restoredefaultpath;
    path(ct_path, path);

    % Add search paths for toolboxes
    if isa(ver, 'char')
        if strcmp(ver, 'legacy')
            ctpath_legacy;
        elseif strcmp(ver, 'new')
            ctpath_new;
        else
            error('Interface version must be "legacy" or "new"');
        end
    else
        error('Invalid input');
    end
end
