function SetMatlabToolboxVersion(ver)
    % Sets the version of the Cantera Matlab Interface
    %
    % :param ver:
    %    integer 0 - legacy 1- new OR
    %    string  'legacy'/'new';
    path('ctpath.m');
    if isa(ver, 'double')
        if ver == 0
            ctpath_legacy;
        elseif ver == 1
            ctpath_new;
        else
            error('Interface version must be 0 (legacy) or 1 (new)');
        end
    elseif isa(ver, 'char')
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