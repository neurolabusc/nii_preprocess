function setenv(name,value)
% an overload replacement for setenv() for static paths
% prevents setting the PATH or LD_LIBRARY_PATH environment variables
% for HPC environments, etc.  David Reddy 20201117
    if ~exist('name','var') || ~exist('value','var')
        warning('usage: setenv(name,value) [stub version for static paths]');
    else
        if (strcmp(name,'PATH') || strcmp(name,'LD_LIBRARY_PATH'))
            warning('Stub setenv() for static paths: NOT exporting %s=%s to the environment\n',name, value);
        else
            builtin('setenv',name,value);
        end
    end
end
